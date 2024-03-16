import unittest
from orcawrapper import orca_input
import numpy as np


class MyTestCase(unittest.TestCase):
    def test_init(self):
        oinp = orca_input("01-inp_xyz/01-xtb2.inp")
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



if __name__ == '__main__':
    unittest.main()
