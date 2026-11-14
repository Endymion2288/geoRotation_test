"""Check material transforms against ACTS's independent row-major decoder rule."""
import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from rotate_material import R, rotate, transform
from compare_material_maps import transform as compare_transform


def native_local(encoded, global_position):
    # ACTS AlgebraJsonConverter.cpp initializes Eigen with
    # rot << rdata[0], rdata[1], ..., rdata[8] (logical rows).
    data = encoded.get('rotation')
    matrix = np.eye(3) if data is None else np.array([
        [data[0], data[1], data[2]],
        [data[3], data[4], data[5]],
        [data[6], data[7], data[8]]])
    translation = np.array(encoded.get('translation') or [0., 0., 0.])
    return np.linalg.solve(matrix, global_position-translation)


class MaterialRotationTest(unittest.TestCase):
    def test_rotated_lookup_preserves_asymmetric_physical_point(self):
        original = {'rotation': None, 'translation': [0., 0., -1822.1075]}
        material = {'binUtility': {'transform': original},
                    'data': [[{'material': [1., 2., 3., 4., 5.], 'thickness': 7.}]]}
        source = {'Surfaces': {'entries': [dict(volume=2, layer=6, approach=1,
                    value={'material': material})]}, 'Volumes': {'entries': []}}
        output = rotate(source)
        value = output['Surfaces']['entries'][0]['value']['material']
        point = np.array([-60.8920685443227, -78.8851533335601, -1820.1075])
        np.testing.assert_allclose(native_local(value['binUtility']['transform'], R@point),
                                   native_local(original, point), atol=1e-12)
        self.assertEqual(value['data'], material['data'])
        self.assertEqual(source['Surfaces']['entries'][0]['value']['material']['binUtility']['transform'], original)

    def test_nonidentity_input_and_comparison_follow_native_decoder(self):
        encoded = {'rotation': [0., -1., 0., 1., 0., 0., 0., 0., 1.],
                   'translation': [3., 5., 7.]}
        point = np.array([11., 13., 17.])
        for decode in [transform, lambda x: compare_transform({'transform': x})]:
            matrix, translation = decode(encoded)
            np.testing.assert_allclose(matrix.T@(point-translation), native_local(encoded, point))

    def test_json_null_translation_is_identity(self):
        for decode in [transform, lambda x: compare_transform({'transform': x})]:
            matrix, translation = decode({'rotation': None, 'translation': None})
            np.testing.assert_array_equal(matrix, np.eye(3))
            np.testing.assert_array_equal(translation, np.zeros(3))


if __name__ == '__main__':
    unittest.main()
