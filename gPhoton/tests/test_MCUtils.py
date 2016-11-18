from __future__ import absolute_import, division, print_function
import unittest
import gPhoton.MCUtils as mc

class TestMCUtilsFunctions(unittest.TestCase):
    def setUp(self):
        self.bands = ['NUV','FUV']
        self.NUV = 'NUV'
        self.FUV = 'FUV'
        self.ra0 = 176.919525856
        self.dec0 = 0.255696872807
        self.t0 = 766525332.995
        self.t1 = 866526576.995
        self.radius = 0.004
        self.eclipse = 23456
        self.maglimit = 30.
        self.detsize = 1.25
        self.xr = [200,400]
        self.yr = [300,500]

    def test_angularSeparation(self):
            self.assertAlmostEqual(
                mc.angularSeparation(10,20,10.1,20.1),0.13720278279273748)


suite = unittest.TestLoader().loadTestsFromTestCase(TestMCUtilsFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)
