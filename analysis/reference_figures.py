#!/usr/bin/env python3
"""Merge source-macro TH1 counts and render native Calypso reference comparisons."""
import argparse
from collections import Counter
import csv
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile

TRUTHS = ('native','production','common_native','common_production')
COMPONENTS = ('smt','flt','prt')
LAYERS = [(s,l) for s in [1,2,3] for l in [0,1,2]]


def histogram_name(component, station, layer):
    return f'h_{component}_{station}_{layer}'


def worker_histogram(file, component, station, layer):
    name=histogram_name(component,station,layer)
    if file.GetKey(name):
        return file.Get(name),False
    # First worker snapshot used TH1::Write() without a key argument. ROOT
    # retained the source object's hist_s_l name, with S/F/P cycles 1/2/3.
    legacy=f'hist_{station}_{layer};{component+1}'
    hist=file.Get(legacy)
    assert hist and hist.GetLineColor()==[600,632,416][component],(file.GetName(),legacy)
    return hist,True


def merge(root, cases):
    import ROOT
    destination = root/'reductions'
    destination.mkdir(exist_ok=True)
    coverage = dict(complete=True, cases={}, source_macro='GitHub ckf_residual_perLayer.C')
    for case in cases:
        tasks = json.loads((root/(case+'_tasks.json')).read_text())
        seen = set()
        counts = Counter()
        per_layer = {f'{s}_{l}': Counter() for s,l in LAYERS}
        maximum = 0.
        examples = []
        histograms = {}
        resources = []
        legacy_files=set()
        for task in tasks:
            identity = (task['input'],task['skip'],task['events'])
            assert identity not in seen
            seen.add(identity)
            source = Path(task['output'])
            assert (source/'exit_code.txt').read_text().strip() == '0', source
            assert json.loads((source/'task.json').read_text()) == task
            observed = json.loads((source/'truth_comparison.json').read_text())
            assert observed['counts']['source_events'] == task['events']
            counts.update(observed['counts'])
            for layer, value in observed['per_layer'].items():
                per_layer[layer].update(value)
            maximum = max(maximum, observed['max_relative_link_vs_production_qop'])
            examples.extend(dict(file_index=task['file_index'], **x) for x in observed['differing_examples'])
            run = json.loads((source/'run.json').read_text())
            assert run['success'] and not any(run[k] for k in ['audit','reset_event_maps','isolate_smoothing'])
            resources.append(dict(file_index=task['file_index'], wall_s=run['wall_s'],
                **json.loads((source/'resources.json').read_text())))
            for truth in TRUTHS:
                file = ROOT.TFile.Open(str(source/(truth+'_3station_counts.root')))
                assert file and not file.IsZombie()
                for component in range(3):
                    for station, layer in LAYERS:
                        name = histogram_name(component,station,layer)
                        hist,legacy = worker_histogram(file,component,station,layer)
                        if legacy:
                            legacy_files.add(file.GetName())
                        assert hist and hist.GetNbinsX()==100
                        assert abs(hist.GetXaxis().GetXmin()+1.2)<1e-12
                        assert abs(hist.GetXaxis().GetXmax()-1.2)<1e-12
                        key = (truth,name)
                        if key not in histograms:
                            clone = hist.Clone(case+'_'+truth+'_'+name)
                            clone.SetDirectory(0)
                            histograms[key] = clone
                        else:
                            histograms[key].Add(hist)
                file.Close()
        for truth in TRUTHS:
            name = case+'_'+truth+'_3station.root'
            with tempfile.TemporaryDirectory(prefix='reference_counts_') as scratch:
                compact = Path(scratch)/name
                output = ROOT.TFile(str(compact),'RECREATE')
                assert not output.IsZombie()
                for (reference,key), hist in histograms.items():
                    if reference == truth:
                        assert hist.Write(key) > 0
                assert not output.TestBit(ROOT.TFile.kWriteError)
                output.Close()
                shutil.copy2(compact,destination/name)
        flow_checks = []
        for truth in TRUTHS:
            for component,stage in enumerate(COMPONENTS):
                for station,layer in LAYERS:
                    value = per_layer[f'{station}_{layer}']
                    prefix = stage+'_'+truth
                    hist = histograms[(truth,histogram_name(component,station,layer))]
                    finite_visible = value[prefix+'_finite']-value[prefix+'_underflow']-value[prefix+'_overflow']
                    # TTreeFormula maps reco q/p=0 or +/-inf to residual=-1.
                    # Those entries remain in the source macro but are excluded
                    # from the separately reported finite-reconstruction metrics.
                    invalid_reco_visible = hist.Integral()-finite_visible
                    assert hist.GetEntries() == value['states']
                    assert 0 <= invalid_reco_visible <= value[prefix+'_invalid_reco']
                    assert hist.GetBinContent(0) == value[prefix+'_underflow']
                    assert hist.GetBinContent(101) == value['states']-value[prefix+'_finite']+value[prefix+'_overflow']-invalid_reco_visible
                    flow_checks.append(dict(truth=truth,stage=stage,station=station,layer=layer,
                        states=value['states'],finite_residual=value[prefix+'_finite'],
                        visible_invalid_reconstruction_entries=int(invalid_reco_visible),
                        entries_and_flow_closure=True))
        comparison = {}
        for first, second in [('native','production'),('common_native','common_production')]:
            difference = {}
            for component in range(3):
                for station, layer in LAYERS:
                    name = histogram_name(component,station,layer)
                    a,b = (histograms[(truth,name)] for truth in [first,second])
                    delta = [b.GetBinContent(i)-a.GetBinContent(i) for i in range(102)]
                    difference[name] = dict(different_bins=sum(x!=0 for x in delta),
                        sum_absolute_bin_difference=sum(abs(x) for x in delta),
                        maximum_absolute_bin_difference=max(abs(x) for x in delta),
                        first_entries=a.GetEntries(), second_entries=b.GetEntries(),
                        first_visible_entries=a.Integral(), second_visible_entries=b.Integral(),
                        first_overflow=a.GetBinContent(101), second_overflow=b.GetBinContent(101))
            comparison[first+'_vs_'+second] = difference
        summary = dict(counts=dict(counts),per_layer={k:dict(v) for k,v in per_layer.items()},
            max_relative_link_vs_production_qop=maximum, differing_examples=examples[:20],
            histogram_comparisons=comparison, resources=resources,
            histogram_flow_checks=flow_checks,
            legacy_worker_key_files=sorted(legacy_files),
            legacy_key_mapping='hist_station_layer cycles 1/2/3 = smt/flt/prt; line colors checked')
        (destination/(case+'_truth_summary.json')).write_text(json.dumps(summary,indent=2)+'\n')
        coverage['cases'][case] = dict(files=len(seen),source_events=counts['source_events'],
            writer_tracks=counts['writer_tracks'], github_selected_tracks=counts['github_selected_tracks'])
        fields=['case','station','layer','stage','truth','states','finite','tail_half','tail_percent',
                'mean_including_overflow','rms_including_overflow','invalid_truth','invalid_reco','overflow']
        with (destination/(case+'_layer_metrics.csv')).open('w') as stream:
            writer = csv.DictWriter(stream,fieldnames=fields)
            writer.writeheader()
            for (station,layer) in LAYERS:
                value=per_layer[f'{station}_{layer}']
                for stage in COMPONENTS:
                    for truth in TRUTHS:
                        key=stage+'_'+truth
                        n=value[key+'_finite']
                        mean=value[key+'_sum']/n if n else None
                        rms=math.sqrt(max(0.,value[key+'_sum2']/n-mean*mean)) if n else None
                        writer.writerow(dict(case=case,station=station,layer=layer,stage=stage,truth=truth,
                            states=value['states'],finite=n,tail_half=value[key+'_below_minus_half'],
                            tail_percent=100*value[key+'_below_minus_half']/n if n else None,
                            mean_including_overflow=mean,rms_including_overflow=rms,
                            invalid_truth=value[key+'_invalid_truth'],invalid_reco=value[key+'_invalid_reco'],
                            overflow=value[key+'_overflow']))
    comparisons={}
    for first,second in [('beforeGeo','afterGeo'),('afterGeo','afterGeo_strict'),
                         ('afterGeo_f00','afterGeo_mask012'),('afterGeo_f00','afterGeo_4station')]:
        if first not in cases or second not in cases:
            continue
        a=read_histograms(destination/(first+'_native_3station.root'),'check_'+first)
        b=read_histograms(destination/(second+'_native_3station.root'),'check_'+second)
        pair={}
        for key,hist in a.items():
            other=b[key]
            delta=[other.GetBinContent(i)-hist.GetBinContent(i) for i in range(102)]
            pair[histogram_name(*key)]=dict(different_bins=sum(d!=0 for d in delta),
                maximum_absolute_bin_difference=max(abs(d) for d in delta),
                sum_absolute_bin_difference=sum(abs(d) for d in delta),
                first_entries=hist.GetEntries(),second_entries=other.GetEntries())
        comparisons[first+'_vs_'+second]=dict(histograms=pair,
            first_scope=coverage['cases'][first],second_scope=coverage['cases'][second])
    (destination/'configuration_comparisons.json').write_text(json.dumps(comparisons,indent=2)+'\n')
    (destination/'coverage.json').write_text(json.dumps(coverage,indent=2)+'\n')
    print(json.dumps(coverage),flush=True)


def read_histograms(path, prefix):
    import ROOT
    file=ROOT.TFile.Open(str(path))
    assert file and not file.IsZombie(),path
    histograms={}
    for component in range(3):
        for station,layer in LAYERS:
            name=histogram_name(component,station,layer)
            hist=file.Get(name).Clone(prefix+'_'+name)
            hist.SetDirectory(0)
            histograms[(component,station,layer)]=hist
    file.Close()
    return histograms


def overlay(first, second, prefix, title, labels):
    """Draw source-produced counts in the reference 3x3 P/F/S layout."""
    import ROOT
    canvas=ROOT.TCanvas(prefix.name,title,1200,800)
    canvas.Divide(3,3)
    for panel in range(1,10):
        pad=canvas.GetPad(panel)
        pad.SetPad(pad.GetXlowNDC(),.96*pad.GetYlowNDC(),
                   pad.GetXlowNDC()+pad.GetWNDC(),.96*(pad.GetYlowNDC()+pad.GetHNDC()))
    retained=[]
    for panel,(station,layer) in enumerate(LAYERS):
        canvas.cd(panel+1)
        maximum=max(hist[(component,station,layer)].GetMaximum()
            for hist in [first,second] for component in range(3))
        leg=ROOT.TLegend(.50,.57,.98,.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(.039)
        for reference,histograms in enumerate([first,second]):
            for component,color in enumerate([ROOT.kBlue,ROOT.kRed,ROOT.kGreen]):
                hist=histograms[(component,station,layer)]
                hist.SetTitle(f'Station {station}, layer {layer}')
                hist.GetXaxis().SetTitle('(P^{reco} - P^{truth})/P^{truth}')
                hist.GetYaxis().SetTitle('Entries')
                hist.SetLineColor(color)
                hist.SetLineStyle(1 if reference==0 else 2)
                hist.SetLineWidth(2)
                hist.SetMinimum(0)
                hist.SetMaximum(1.2*maximum)
                hist.Draw('hist' if reference==component==0 else 'hist same')
                leg.AddEntry(hist,labels[reference]+' '+COMPONENTS[component],'l')
        leg.Draw()
        latex=ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(.05)
        latex.DrawLatex(.15,.85,f'S{station} L{layer}')
        retained.extend([leg,latex])
    canvas.cd(0)
    latex=ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(.018)
    latex.SetTextAlign(23)
    latex.DrawLatex(.5,.988,title)
    canvas.SaveAs(str(prefix)+'.png')
    canvas.SaveAs(str(prefix)+'.pdf')
    saved=ROOT.TFile(str(prefix)+'.root','RECREATE')
    canvas.Write('reference_canvas')
    saved.Close()
    return canvas,retained+[latex]


def plot(root,cases):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    destination=root/'figures'
    destination.mkdir(exist_ok=True)
    reductions=root/'reductions'
    coverage=json.loads((reductions/'coverage.json').read_text())
    ROOT.gROOT.ProcessLine('.L '+str(root/'reference/adapted/ckf_residual_perLayer.C'))
    data={(case,truth):read_histograms(reductions/(case+'_'+truth+'_3station.root'),case+'_'+truth)
        for case in cases for truth in TRUTHS}
    maxima=[max(hist[(component,station,layer)].GetMaximum()
        for hist in data.values() for component in range(3)) for station,layer in LAYERS]
    for panel,maximum in enumerate(maxima):
        os.environ[f'REFERENCE_YMAX_{panel}']=str(maximum)
    index=[]
    for case in cases:
        scope=coverage['cases'][case]
        for truth in ['native','production']:
            reference='SiHit-linked HepMC truth' if truth=='native' else 'Production truth (same majority barcode)'
            prefix=destination/('110051_'+case+'_'+truth+'_nine_layers')
            os.environ['REFERENCE_SCOPE_LABEL']=f"110051 | {case} | {scope['files']} files | {scope['source_events']:,} source events | {reference}"
            ROOT.ckf_residual_perLayer(str(reductions/(case+'_'+truth+'_3station.root')),
                str(prefix),case,'t_eQOP',True)
            index.append(dict(name=prefix.name,description=os.environ['REFERENCE_SCOPE_LABEL'],reference='geometryRotation.pdf p17; original GitHub macro'))
        title=f"110051 | {case} | {scope['files']} files | {scope['source_events']:,} source events | Only truth reference changes"
        prefix=destination/('110051_'+case+'_truth_overlay')
        canvas,retained=overlay(data[(case,'native')],data[(case,'production')],prefix,title,
            ['SiHit link','Production'])
        index.append(dict(name=prefix.name,description=title,reference='GitHub counts; solid=link, dashed=production'))
        canvas.Close()
        prefix=destination/('110051_'+case+'_common_truth_overlay')
        canvas,retained=overlay(data[(case,'common_native')],data[(case,'common_production')],prefix,
            title+' | Common finite truth states',['SiHit link','Production'])
        index.append(dict(name=prefix.name,description=title+' | Common finite truth states',reference='GitHub counts; common truth coverage'))
        canvas.Close()
    for first,second,name in [('beforeGeo','afterGeo','geometry_overlay'),
                              ('afterGeo','afterGeo_strict','material_overlay')]:
        if first not in cases or second not in cases:
            continue
        scope=coverage['cases'][first]
        title=f"110051 | {scope['files']} files | {scope['source_events']:,} source events | SiHit-linked HepMC truth | {first} vs {second}"
        prefix=destination/('110051_'+name)
        canvas,retained=overlay(data[(first,'native')],data[(second,'native')],prefix,title,[first,second])
        index.append(dict(name=prefix.name,description=title,reference='GitHub-macro histograms; solid=first, dashed=second'))
        canvas.Close()
    # Preserve the individual ROOT vector exports, including their page boxes.
    sys.path.insert(0,str(root/'software_tools'))
    from pypdf import PdfReader,PdfWriter
    writer=PdfWriter()
    for item in index:
        page=PdfReader(destination/(item['name']+'.pdf'))
        assert len(page.pages)==1
        writer.add_page(page.pages[0])
    with (destination/'reference_comparison.pdf').open('wb') as stream:
        writer.write(stream)
    combined = PdfReader(destination/'reference_comparison.pdf')
    assert len(combined.pages) == len(index)
    pages = []
    for number,item in enumerate(index):
        original = PdfReader(destination/(item['name']+'.pdf')).pages[0]
        page = combined.pages[number]
        assert page.get_contents().get_data() == original.get_contents().get_data()
        assert list(page.cropbox) == list(original.cropbox)
        assert list(page.mediabox) == list(original.mediabox)
        assert page.rotation == original.rotation
        pages.append(dict(name=item['name'],content_and_page_boxes_preserved=True))
    (destination/'figure_validation.json').write_text(json.dumps(dict(pages=pages,
        page_count=len(pages),pdf_assembler='pypdf'),indent=2)+'\n')
    (destination/'figure_index.json').write_text(json.dumps(index,indent=2)+'\n')
    (destination/'figure_index.md').write_text('# 110051 reference comparisons\n\n'+
        '\n'.join(f"- [{item['name']}]({item['name']}.png): {item['description']}" for item in index)+'\n')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--phase',choices=['merge','plot','all'],default='all')
    parser.add_argument('--cases',nargs='+',default=['beforeGeo','afterGeo','afterGeo_strict'])
    parser.add_argument('--plot-cases',nargs='+',help='Subset to plot after merging the complete control matrix')
    args=parser.parse_args()
    if args.phase in ['merge','all']:
        merge(args.root,args.cases)
    if args.phase in ['plot','all']:
        plot(args.root,args.plot_cases or args.cases)


if __name__=='__main__':
    main()
