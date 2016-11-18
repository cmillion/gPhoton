from __future__ import absolute_import, division, print_function
import unittest
import gPhoton.dbasetools as dbt

class TestDbasetoolsFunctions(unittest.TestCase):

    def setUp(self):
        self.bands = ['NUV','FUV']

    def test_fGetTimeRanges(self):
        self.assertTrue(True)

suite = unittest.TestLoader().loadTestsFromTestCase(TestDbasetoolsFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)
