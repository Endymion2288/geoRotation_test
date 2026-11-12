#!/usr/bin/env python3
"""Adapt only I/O, labels and unused crashing display blocks in the GitHub macro."""
import argparse
import difflib
from pathlib import Path
import re


def adapt(source):
    text = source
    begin = text.index('void ckf_residual_perLayer(')
    end = text.index(') {', begin)+3
    text = text[:begin]+'''void ckf_residual_perLayer(const char* inputPath, const char* outputPrefix,
             const char* label, const char* truthBranch="t_eQOP", bool prebinned=false,
             int station=1, int layer=0, int drawType=1,
             bool doFit=false, double minY=0, double maxY=1) {'''+text[end:]
    text, n = re.subn(r'  std::vector<std::string>files = \{.*?\n  \};',
        '  std::vector<std::string>files = {inputPath};', text, count=1, flags=re.S)
    assert n == 1
    text = text.replace('abs(t_eQOP)', 'abs(reference_truth)')
    line = '    std::cout<<"events  = "<<tree->GetEntries()<<std::endl;'
    assert line in text
    text = text.replace(line, '''    if (tree) {
      tree->AddFriend("production", "production.root");
      tree->SetAlias("reference_truth", truthBranch);
      std::cout << "tracks = " << tree->GetEntries() << std::endl;
    }''')
    draw = '''        tree->Draw(Form("%s>>hist_%i_%i", format.second.c_str(), i, j),
                   Form("volume_id==%i&&layer_id==%i&&chi2Sum/"
                        "NDF<20&&rMax_smt<90&&NDF>9",
                        station, layer));'''
    assert draw in text
    text = text.replace(draw, '''        if (prebinned) {
          int component = format.second.find("smt") != std::string::npos ? 0 :
                          format.second.find("flt") != std::string::npos ? 1 : 2;
          auto saved = static_cast<TH1*>(tfile->Get(Form("h_%d_%d_%d", component, i, j)));
          if (!saved) throw std::runtime_error("Missing reference histogram");
          hist->Add(saved);
        } else {
'''+draw+'''
        }''')
    # The active source has num=3 but later accesses hists[3] and numStations[2].
    # Its completed nine-panel cCompF is the figure on PDF page 17.
    text = text[:text.index('    TCanvas* cCompB =')]
    text = text.replace('"beforeGeo_110044(μ+)_smt"', 'Form("%s smt",label)')
    text = text.replace('"beforeGeo_110044(μ+)_flt"', 'Form("%s flt",label)')
    text = text.replace('"beforeGeo_110044(μ+)_prt"', 'Form("%s prt",label)')
    text = text.replace('  cCompF->Divide(3,3);', '''  cCompF->Divide(3,3);
  for (int panel=1;panel<=9;++panel) {
    auto pad = static_cast<TPad*>(cCompF->GetPad(panel));
    pad->SetPad(pad->GetXlowNDC(), .96*pad->GetYlowNDC(),
                pad->GetXlowNDC()+pad->GetWNDC(), .96*(pad->GetYlowNDC()+pad->GetHNDC()));
  }''')
    text = text.replace('      h1->SetMaximum(1.2*ymax);', '''      const char* sharedMaximum = gSystem->Getenv(Form("REFERENCE_YMAX_%d",ilay));
      if (sharedMaximum) ymax = std::atof(sharedMaximum);
      h1->SetMaximum(1.2*ymax);''')
    text += '''
  cCompF->cd(0);
  TLatex scope;
  scope.SetNDC(); scope.SetTextSize(0.018); scope.SetTextAlign(23);
  const char* scopeLabel = gSystem->Getenv("REFERENCE_SCOPE_LABEL");
  scope.DrawLatex(0.5,0.988,scopeLabel ? scopeLabel : label);
  cCompF->SaveAs((std::string(outputPrefix)+".png").c_str());
  cCompF->SaveAs((std::string(outputPrefix)+".pdf").c_str());
  TFile output((std::string(outputPrefix)+".root").c_str(),"RECREATE");
  for (int component=0;component<3;++component)
    for (int s=0;s<3;++s) for(int l=0;l<3;++l)
      hists[component][s][l]->Write(Form("h_%d_%d_%d",component,s+1,l));
  cCompF->Write("reference_canvas");
  output.Close();
}
'''
    assert draw in text
    return text


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    original = args.source.read_text()
    adapted = adapt(original)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(adapted)
    args.output.with_suffix('.patch').write_text(''.join(difflib.unified_diff(
        original.splitlines(keepends=True), adapted.splitlines(keepends=True),
        fromfile='GitHub/ckf_residual_perLayer.C', tofile='adapted/ckf_residual_perLayer.C')))


if __name__ == '__main__':
    main()
