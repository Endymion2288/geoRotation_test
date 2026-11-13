"""A beforeGeo worker must not supply geometry to an afterGeo subprocess."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from run_slope_simulation import reconstruction_environment


class SlopeEnvironmentTests(unittest.TestCase):
    def test_cross_variant_dependency_fallback(self):
        libraries = dict(beforeGeo='/tmp/chart/software/lib', afterGeo='/tmp/after/software/lib')
        inherited = dict(GEOROTATION_LIBRARY_DIR=libraries['beforeGeo'],
                         LD_LIBRARY_PATH='/tmp/chart/software/lib:/cvmfs/calypso/lib:/cvmfs/athena/lib',
                         PYTHONPATH='/tmp/chart/software/python:/cvmfs/calypso/python')
        env = reconstruction_environment(Path(libraries['afterGeo']), libraries, inherited)
        self.assertEqual(env['LD_LIBRARY_PATH'].split(':'),
                         ['/tmp/after/software/lib', '/cvmfs/calypso/lib', '/cvmfs/athena/lib'])
        self.assertEqual(env['PYTHONPATH'].split(':'),
                         ['/tmp/after/software/python', '/cvmfs/calypso/python'])
        self.assertEqual(inherited['GEOROTATION_LIBRARY_DIR'], libraries['beforeGeo'])

    def test_remove_parent_overlay_even_when_not_a_target_mode(self):
        inherited = dict(GEOROTATION_LIBRARY_DIR='/tmp/old/software/lib',
                         LD_LIBRARY_PATH='/tmp/old/software/lib:/tmp/after/software/lib:/cvmfs/base/lib',
                         PYTHONPATH='/tmp/old/software/python:/cvmfs/base/python')
        env = reconstruction_environment(Path('/tmp/after/software/lib'),
                                         dict(afterGeo='/tmp/after/software/lib'), inherited)
        self.assertEqual(env['LD_LIBRARY_PATH'], '/tmp/after/software/lib:/cvmfs/base/lib')
        self.assertEqual(env['PYTHONPATH'], '/tmp/after/software/python:/cvmfs/base/python')


if __name__ == '__main__':
    unittest.main()
