from __future__ import absolute_import, division, print_function
import unittest
import gPhoton.dbasetools as dbt

class TestDbasetoolsFunctions(unittest.TestCase):

    def setUp(self):
        self.bands = ['NUV','FUV']

    def test_compute_exptime_1(self):
        skypos = (170.04661625965298, 9.85073479305886)
        trange = (891394895.0, 891394925.0)
        self.assertAlmostEqual(dbt.compute_exptime('NUV',trange),
                               dbt.compute_exptime('NUV',trange,skypos=skypos),
                               places=12)

    def test_compute_exptime_2(self):
        skypos = [176.919525856024,0.255696872807351]
        trange = [766525332.995,766526576.995]
        self.assertAlmostEqual(dbt.compute_exptime('NUV',trange),
                               dbt.compute_exptime('NUV',trange,skypos=skypos),
                               places=12)

    def test_compute_exptime_3(self):
        """Check that exposure time over a time range and exposure time
        over a time range passed through fGetTimeRanges return the same."""

suite = unittest.TestLoader().loadTestsFromTestCase(TestDbasetoolsFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)
