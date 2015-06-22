import unittest
import gPhoton.gMap as gm
import gPhoton.gFind as gf
import gPhoton.gAperture as ga

"""Regression tests for the command line utilities."""
class TestRegression(unittest.TestCase):
    def setUp(self):
        self.bands = ['NUV','FUV']
        self.skypos = [176.919525856024,0.255696872807351]
        self.tranges = [[ 765579107.995, 765579821.995 ],
                        [ 766525332.995, 766526576.995 ],
                        [ 919754986.995, 919755916.995 ]]
        self.trange = self.tranges[0]
        self.radius = 0.01
        self.annulus = [0.02,0.03]
        self.parser = ga.setup_parser()
        self.args = self.parser.parse_args()

    def test_blah(self):
        self.assertTrue(True)

suite = unittest.TestLoader().loadTestsFromTestCase(TestRegression)
unittest.TextTestRunner(verbosity=2).run(suite)
