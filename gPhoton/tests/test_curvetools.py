from __future__ import absolute_import, division, print_function
import unittest
import gPhoton.curvetools as ct

class TestCurvetoolsFunctions(unittest.TestCase):

    def setUp(self):
        self.bands = ['NUV','FUV']

    # def test_out_of_bounds_event_1(self):
    #     # Test an unusual edge case bug where a photon event might be returned
    #     #  that is outside of the requested time range
    #     events = ct.query_photons('NUV', 179.85709, -21.97462,
    #         [[823515332.99998, 823515445.9950199]], 0.025)
    #     self.assertTrue(events['t'].min()>=823515332.99998)

    def test_out_of_bounds_event_2(self):
        data = ct.get_curve('NUV', 179.85709, -21.97462, 0.025,
            trange=[823515332.99998, 823515445.9950199], detsize=1.25)
        self.assertTrue(data['photons']['t'].min()>=data['params']['trange'][0])

suite = unittest.TestLoader().loadTestsFromTestCase(TestCurvetoolsFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)
