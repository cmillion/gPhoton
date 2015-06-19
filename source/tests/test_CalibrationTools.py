import CalibrationTools as CalTools
import unittest
import numpy as np
import os

class TestCalibrationToolsFunctions(unittest.TestCase):

    def setUp(self):
        self.FUV = 'FUV'
        self.NUV = 'NUV'

    def test_load_txy(self):
        # You have to make a fake csv file or else do a full gPhoton
        self.assertTrue(1==1)

    def test_compute_deadtime(self):
        self.assertTrue(1==1)

    def test_compute_shutter(self):
        self.assertTrue(1==1)

    def test_compute_exposure(self):
        self.assertTrue(1==1)

    def test_create_rr(self):
        self.assertTrue(1==1)

    def test_write_rr(self):
        self.assertTrue(1==1)

    def test_write_rrhr(self):
        self.assertTrue(1==1)

    def test_write_int(self):
        self.assertTrue(1==1)

#       def tearDown():
#

suite = unittest.TestLoader().loadTestsFromTestCase(TestCalibrationToolsFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)

