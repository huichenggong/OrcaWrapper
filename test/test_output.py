import unittest
from orcawrapper import orca_opt_output
import numpy as np


class MyTestCase(unittest.TestCase):
    def test_check_terminated_normally(self):
        print("# Test : check_terminated_normally(), check if the job terminated normally.")
        oout = orca_opt_output("02-OPT/opt.out")
        self.assertTrue(oout.terminated_normally)

    def test_get_opt_traj(self):
        print("# Test : get_opt_traj(), read optimized trajectory.")
        oout = orca_opt_output("02-OPT/opt.out")
        np.testing.assert_allclose(oout.xyz[0], np.array([[ 0.000000,   0.000000,   0.000000],
                                                          [ 1.250000,   0.000000,   0.000000],
                                                          [-0.587148,   0.939049,   0.000000],
                                                          [-0.587148,  -0.939049,  -0.000000],
                                                          ])
                                   )
        np.testing.assert_allclose(oout.xyz[1], np.array([[ 0.011101,  -0.000000,  -0.000000],
                                                          [ 1.220982,   0.000000,   0.000000],
                                                          [-0.578190,   0.946327,   0.000000],
                                                          [-0.578190,  -0.946327,   0.000000],
                                                          ])
                                   )
        self.assertEqual(len(oout.xyz), 5)
        self.assertAlmostEqual(oout.time, 30.620)
        self.assertDictEqual(oout.basic_info, {'n_atoms': 4, 'n_basis_functions': 38})

    def test_get_opt_info(self):
        print("# Test : get_opt_info(), read optimized information.")
        oout = orca_opt_output("02-OPT/opt.out")
        print(oout.get_opt_info())


if __name__ == '__main__':
    unittest.main()
