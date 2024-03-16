import unittest
from orcawrapper import orca_input
import numpy as np


class MyTestCase(unittest.TestCase):
    def test_init(self):
        print("# Test : __init__(), read correct coordinate.")
        oinp = orca_input("01-inp_xyz/01-xtb2.inp")
        print(oinp)
        np.testing.assert_allclose(oinp.xyz[:6,:],
                                   np.array([[0.824, 8.697, 1.039],
                                             [0.561, 9.629, 1.288],
                                             [1.260, 8.253, 1.822],
                                             [ 1.465, 8.725, 0.272],
                                             [-0.417, 7.910, 0.635],
                                             [-0.332, 6.725, 0.286],
                                             ]))
        np.testing.assert_allclose(oinp.xyz[-2:, :],
                                   np.array([[0.87810201,    0.84350579,  -6.866131],
                                             [1.40472737, 1.99094572, -5.93680473]]))
        self.assertEqual(oinp.xyz.shape, (500, 3))

    def test_update_xyz(self):
        print("# Test : update_xyz(), update the coordinate.")
        oinp = orca_input("01-inp_xyz/01-xtb2.inp")
        new_xyz = np.concatenate(([[1, 2, 3],
                                   [4, 5, 6],
                                   [7, 8, 9],
                                   ], np.random.rand(497, 3)))
        oinp.update_xyz(new_xyz)
        np.testing.assert_allclose(oinp.xyz, new_xyz)
        np.testing.assert_allclose(oinp.xyz[:3, :],
                                   np.array([[1, 2, 3],
                                             [4, 5, 6],
                                             [7, 8, 9]]))

    def test_write(self):
        print("# Test : write(), write an input file.")
        oinp = orca_input("01-inp_xyz/01-xtb2.inp")
        # set random seed
        np.random.seed(0)
        new_xyz = np.concatenate(([[1, 2, 3],
                                   [4, 5, 6],
                                   [7, 8, 9],
                                   ], np.random.rand(497, 3)))
        oinp.update_xyz(new_xyz)
        oinp.write("01-inp_xyz/01-xtb2_new.inp")
        oinp_new = orca_input("01-inp_xyz/01-xtb2_new.inp")
        np.testing.assert_allclose(oinp_new.xyz, new_xyz)
        np.testing.assert_allclose(oinp_new.xyz[:3, :],
                                   np.array([[1, 2, 3],
                                             [4, 5, 6],
                                             [7, 8, 9]]))



if __name__ == '__main__':
    unittest.main()
