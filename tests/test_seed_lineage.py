"""ROOT regression for duplicate event IDs and same-barcode competing seeds."""
from array import array
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from native_purity import attach_final_purity


@unittest.skipUnless(importlib.util.find_spec('ROOT'), 'PyROOT required')
class SeedLineageTest(unittest.TestCase):
    def test_actual_seed_and_source_entry_control_the_gate(self):
        import ROOT
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = ROOT.TFile(str(root/'track_states.root'), 'RECREATE')
            tree = ROOT.TTree('tree', 'native writer fixture')
            scalars = {}
            for name, value in dict(event_nr=1, track_nr=0, seed_nr=0,
                                    majorityParticleId=13, NDF=12).items():
                scalars[name] = array('i', [value])
                tree.Branch(name, scalars[name], name+'/I')
            for name, value in dict(chi2Sum=1., rMax_smt=10.).items():
                scalars[name] = array('d', [value])
                tree.Branch(name, scalars[name], name+'/D')
            vectors = {}
            for name, value in dict(t_eQOP=.00001, volume_id=25., layer_id=2.,
                    eQOP_prt=.01, eQOP_flt=.01, eQOP_smt=.01).items():
                vectors[name] = ROOT.std.vector('double')()
                vectors[name].push_back(value)
                tree.Branch(name, vectors[name])
            tree.Fill()
            tree.Fill()
            tree.Write()
            output.Close()

            seed_path = root/'seeds.root'
            output = ROOT.TFile(str(seed_path), 'RECREATE')
            tree = ROOT.TTree('tree', 'seed writer fixture')
            event = array('i', [1])
            tree.Branch('event_number', event, 'event_number/I')
            vectors = {name: ROOT.std.vector('double')() for name in
                ['seedIndex', 'majorityParticleId', 'track_QoverP', 't_QoverP',
                 'track_theta', 't_theta']}
            for name, vector in vectors.items():
                tree.Branch(name, vector)
            # Both events have the same production ID and both seeds the same
            # barcode. Seed 0 fails in entry 0, passes in entry 1; seed 1 does
            # the reverse. Only exact source-entry + seed association is valid.
            for momenta in [[.1, .01], [.01, .1]]:
                values = dict(seedIndex=[0, 1], majorityParticleId=[13, 13],
                    track_QoverP=momenta, t_QoverP=[.01, .01],
                    track_theta=[.2, .2], t_theta=[.2, .2])
                for name, vector in vectors.items():
                    vector.clear()
                    for value in values[name]:
                        vector.push_back(value)
                tree.Fill()
            tree.Write()
            output.Close()
            rows = [dict(source_entry=i, event=1,
                final=[dict(clusters=18, majority_hits=18, majority_barcode=13)])
                for i in range(2)]
            (root/'production_truth.jsonl').write_text(
                ''.join(json.dumps(row)+'\n' for row in rows))
            summary = attach_final_purity(root, seed_path)
            self.assertEqual(summary['counts']['baseline_selected_tracks'], 2)
            self.assertEqual(summary['counts']['selected_tracks'], 1)
            self.assertEqual(summary['seed_counts']['barcode_groups_with_multiple_seeds'], 2)
            # Also catches gDirectory changes caused by opening the seed ROOT.
            source = ROOT.TFile.Open(str(root/'production.root'))
            friend = source.Get('production')
            self.assertTrue(friend)
            self.assertEqual([int(row.seed_good) for row in friend], [0, 1])
            source.Close()


if __name__ == '__main__':
    unittest.main()
