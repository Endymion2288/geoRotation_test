#!/usr/bin/env python3
"""Check actual CKF2 default/explicit ON/OFF output on a frozen small sample."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO/'analysis'))


def snapshot(path):
    import ROOT
    scalar = ['event_nr', 'track_nr', 'nStates', 'nMeasurements', 'nHoles',
              'nOutliers', 'chi2Sum', 'NDF', 'majorityParticleId']
    vectors = ['volume_id', 'layer_id', 'module_id', 'index_hit', 't_eQOP',
               'l_x_hit', 'l_y_hit', 'g_x_hit', 'g_y_hit', 'g_z_hit', 'pathLength']
    vectors += [prefix+parameter+'_'+stage for prefix in ['', 'err_']
                for parameter in ['eLOC0', 'eLOC1', 'ePHI', 'eTHETA', 'eQOP', 'eT']
                for stage in ['prt', 'flt', 'smt']]
    columns = ROOT.RDataFrame('tree', str(path)).AsNumpy(scalar+vectors)
    result = {name: columns[name] for name in scalar}
    for name in vectors:
        values = [np.asarray(row) for row in columns[name]]
        result[name] = np.concatenate(values)
        result[name+'_offsets'] = np.cumsum([0]+[len(row) for row in values])
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference-task', type=Path, required=True)
    parser.add_argument('--software', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--previous-on', type=Path, required=True)
    parser.add_argument('--reuse-scratch', type=Path,
        help='Resume reduction of already successful reconstruction after a plotting failure')
    args = parser.parse_args()
    task = json.loads(args.reference_task.read_text())[0]
    assert task['events'] <= 20 and task['variant'] == 'afterGeo'
    root = args.output.resolve()
    root.mkdir(exist_ok=args.reuse_scratch is not None, parents=True)
    (root/'input_task.json').write_text(json.dumps(task, indent=2)+'\n')
    library = args.software.resolve()/'lib'
    assert os.environ['GEOROTATION_LIBRARY_DIR'] == str(library)
    modes = {'default': [], 'explicit_on': ['--isolate-smoothing'],
             'legacy_off': ['--shared-smoothing']}
    # Retain worker scratch on a reduction failure so tracking is not repeated.
    with nullcontext(args.reuse_scratch or Path(tempfile.mkdtemp(
            prefix='georotation_smoothing_smoke_'))) as temporary:
        scratch = Path(temporary)
        (root/'scratch_path.txt').write_text(str(scratch)+'\n')

        def reconstruct(item):
            name, flags = item
            directory = scratch/name
            destination = root/name
            if args.reuse_scratch:
                run = json.loads((destination/'run.json').read_text())
                assert run['success'] and (directory/'track_states.root').is_file()
                return name
            directory.mkdir()
            destination.mkdir()
            command = ['python3', str(REPO/'configs/reference_tracking.py'),
                '--input', task['input'], '--events', str(task['events']),
                '--skip', str(task['skip']), '--material', task['material'],
                '--variant', task['variant'], '--mask', task['mask'],
                '--native-final-purity', *flags,
                '--initial-variance-inflation', *map(str, task['initial_variance_inflation'])]
            (destination/'command.json').write_text(json.dumps(command, indent=2)+'\n')
            with (directory/'reconstruction.log').open('w') as log:
                status = subprocess.run(['/usr/bin/time', '-v', '-o', str(destination/'resource.txt'),
                    *command], cwd=directory, stdout=log, stderr=subprocess.STDOUT).returncode
            with (directory/'reconstruction.log').open('rb') as source, gzip.open(
                    destination/'reconstruction.log.gz', 'wb') as target:
                shutil.copyfileobj(source, target)
            if (directory/'run.json').exists():
                shutil.copy2(directory/'run.json', destination/'run.json')
            assert status == 0, (name, status)
            run = json.loads((destination/'run.json').read_text())
            assert run['success'] and run['isolate_smoothing'] == (name != 'legacy_off')
            assert {str(library/('lib'+target+'.so')) for target in [
                'FaserActsKalmanFilter', 'FaserActsGeometry', 'FaserActsGeometryLib']
                } <= set(run['loaded_acts_libraries'])
            return name

        # Two independent single-threaded Athena processes, then ROOT reduction.
        with ThreadPoolExecutor(max_workers=2) as pool:
            list(pool.map(reconstruct, modes.items()))
        import ROOT
        ROOT.gROOT.SetBatch(True)
        ROOT.gROOT.ProcessLine('.L '+task['macro'])
        data = {}
        histograms = {}
        from reference_figures import read_histograms, overlay
        from native_purity import attach_final_purity
        for name in modes:
            directory = scratch/name
            os.chdir(directory)
            data[name] = snapshot(directory/'track_states.root')
            np.savez_compressed(root/name/'writer_snapshot.npz', **data[name])
            (directory/'reference_3station_forward.root').unlink(missing_ok=True)
            statistics = attach_final_purity(directory)
            assert statistics['counts']['source_events'] == task['events']
            assert statistics['counts']['writer_tracks'] == len(data[name]['track_nr'])
            os.environ['REFERENCE_SCOPE_LABEL'] = (
                f"110051 | afterGeo | 1 file | entries [{task['skip']},{task['skip']+task['events']}) | {name}")
            ROOT.ckf_residual_perLayer(str(directory/'reference_3station_forward.root'),
                str(directory/'native'), 'afterGeo', 't_eQOP', False)
            source = ROOT.TFile.Open(str(directory/'native.root'))
            output = ROOT.TFile(str(root/name/'native_counts.root'), 'RECREATE')
            for component in range(3):
                for station in [1, 2, 3]:
                    for layer in [0, 1, 2]:
                        hist = source.Get(f'h_{component}_{station}_{layer}')
                        assert hist and hist.Write(f'h_{component}_{station}_{layer}') > 0
            output.Close()
            source.Close()
            histograms[name] = read_histograms(root/name/'native_counts.root', name)
        default = data['default']
        assert default.keys() == data['explicit_on'].keys()
        for key, values in default.items():
            assert np.array_equal(values, data['explicit_on'][key], equal_nan=True), key
            if '_smt' not in key:
                assert np.array_equal(values, data['legacy_off'][key], equal_nan=True), key
        changed = int(np.sum(~np.isclose(default['eQOP_smt'],
            data['legacy_off']['eQOP_smt'], atol=0, rtol=0, equal_nan=True)))
        assert changed > 0, 'Sample must demonstrate the shared-state overwrite'
        previous = read_histograms(args.previous_on/'native_3station_counts.root', 'previous')
        for key, hist in histograms['default'].items():
            for other in [histograms['explicit_on'][key], previous[key]]:
                assert hist.GetEntries() == other.GetEntries(), key
                for index in range(102):
                    assert hist.GetBinContent(index) == other.GetBinContent(index), key
                    assert hist.GetBinError(index) == other.GetBinError(index), key
        figures = root/'figures'
        figures.mkdir()
        counts = {}
        all_counts = {name: {str(key): dict(
            edges=[hist.GetBinLowEdge(i) for i in range(1, 102)],
            sumw=[hist.GetBinContent(i) for i in range(102)],
            sumw2=[hist.GetBinError(i)**2 for i in range(102)],
            raw_entries=hist.GetEntries()) for key, hist in group.items()}
            for name, group in histograms.items()}
        (root/'all_histogram_counts.json').write_text(json.dumps(all_counts, indent=2)+'\n')
        from smoothing_recheck import draw_smoothed
        from pypdf import PdfReader, PdfWriter
        pages = []
        index = []
        title = (f"110051 | afterGeo | 1 file | {task['events']} source entries "
                 f"[{task['skip']}, {task['skip']+task['events']})\n"
                 'Rotated ACTS sensor states | native SiHit-linked HepMC truth | GitHub cuts')
        for other, label in [('explicit_on', 'Explicit isolated ON'), ('legacy_off', 'Legacy OFF')]:
            prefix = figures/('default_vs_'+other)
            counts[other] = draw_smoothed([histograms[other], histograms['default']],
                prefix, title, [label, 'New default'])
            pages.append(prefix.with_suffix('.pdf'))
            index.append(dict(png=str(prefix.with_suffix('.png')), pdf=str(pages[-1])))
        prefix = figures/'default_vs_legacy_PFS'
        overlay(histograms['legacy_off'], histograms['default'], prefix,
            '#splitline{'+title.split('\n')[0]+'}{'+title.split('\n')[1]+'}',
            ['Legacy OFF', 'New default'])
        pages.append(prefix.with_suffix('.pdf'))
        index.append(dict(png=str(prefix.with_suffix('.png')), pdf=str(pages[-1])))
        writer = PdfWriter()
        for page in pages:
            writer.append(PdfReader(str(page)))
        with (figures/'smoothing_default_validation.pdf').open('wb') as output:
            writer.write(output)
        (root/'bin_counts.json').write_text(json.dumps(counts, indent=2)+'\n')
        (root/'figure_index.json').write_text(json.dumps(index, indent=2)+'\n')
        result = dict(source_events=task['events'], input=task['input'], skip=task['skip'],
            writer_tracks=len(default['track_nr']), measurement_states=len(default['eQOP_smt']),
            snapshot_arrays_equal_default_vs_on=len(default),
            unchanged_non_smoothed_arrays_vs_off=True, changed_smoothed_qop=changed,
            previous_on=str(args.previous_on), previous_on_native_histograms_equal=27,
            histogram_bins_including_flow=102,
            raw_entries_and_sumw2_equal=True, software=str(args.software),
            histogram_cuts='GitHub: chi2Sum/NDF<20 && NDF>9 && rMax_smt<90; no purity gate for code regression')
        (root/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
        shutil.rmtree(scratch)
        print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
