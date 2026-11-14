import sys
from pathlib import Path
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from charge_material_noise import FIELDS, ROTATE_F_TO_A, noise, noise_difference, trace_gate, trace_rows


def recorded(variant, rounded_theta=True, scattering=True):
    direction_f = np.array([.003, -.002, np.sqrt(1.-.003**2-.002**2)])
    native = direction_f if variant == 'beforeGeo' else ROTATE_F_TO_A @ direction_f
    theta0 = np.float32(.000117)
    square = float(theta0)**2
    phi = square/(native[0]**2+native[1]**2) if scattering else 0.
    theta = float(np.float32(theta0*theta0)) if rounded_theta else square
    row = dict(zip(FIELDS, (phi, theta if scattering else 0., 0.)))
    row.update({'dirF'+axis:float(value) for axis, value in zip('xyz', direction_f)})
    return row


class NativeNoiseCoordinates(unittest.TestCase):
    def test_disabled_trace_is_unmeasured_and_required_trace_cannot_be_missing(self):
        with tempfile.TemporaryDirectory() as directory:
            old, new = [Path(directory)/name for name in ('old', 'new')]
            old.mkdir()
            new.mkdir()
            self.assertFalse(trace_gate(old, new, False)['trace_recorded'])
            self.assertIsNone(trace_rows(new))
            with self.assertRaises(FileNotFoundError):
                trace_gate(old, new, True)
            (new/'material_trace.csv').write_text('source_entry,candidate,target\n')
            with self.assertRaises(AssertionError):
                trace_gate(old, new, False)

    def test_isotropic_noise_is_invariant_under_the_full_frame_rotation(self):
        a, b = [noise(recorded(v, rounded_theta=False), v) for v in ('beforeGeo', 'afterGeo')]
        self.assertLess(noise_difference(a, b)['physical_direction_Q_scaled_max'], 2e-15)

    def test_single_precision_square_is_detected_and_its_tensor_orientation_changes(self):
        a, b = [noise(recorded(v), v) for v in ('beforeGeo', 'afterGeo')]
        for row in (a, b):
            self.assertEqual(row['float_square_relation_error'], 0.)
            self.assertLess(abs(row['double_square_relation_error']), 1e-15)
            self.assertAlmostEqual(row['isotropy_error'], row['expected_isotropy_from_float_square'], places=14)
        delta = noise_difference(a, b)
        self.assertGreater(delta['physical_direction_Q_scaled_max'], 1e-9)
        self.assertLess(delta['same_double_square_Q_scaled_max'], 2e-15)

    def test_zero_material_scattering_remains_a_valid_record(self):
        row = noise(recorded('beforeGeo', scattering=False), 'beforeGeo')
        self.assertFalse(row['scattering'])
        self.assertIsNone(row['isotropy_error'])
        self.assertEqual(noise_difference(row, row)['physical_direction_Q_scaled_max'], 0.)


if __name__ == '__main__':
    unittest.main()
