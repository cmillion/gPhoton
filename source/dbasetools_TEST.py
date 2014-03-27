import dbasetools as dt
import unittest
import numpy as np
import os

# It would be really great to come up with a way to test all of these
#  without assuming that the database is static. If a correction is
#  made to the database at some point, then almost every one of these
#  tests will fail and need to be rewritten.

class TestDbasetoolsFunctions(unittest.TestCase):

	def setUp(self):
		self.foo = 1

	def test_fGetTimeRanges(self):
		self.assertTrue(True)

	def test_compute_exptime(self):
		self.assertTrue(True)

	def test_mcat_skybg(self):
		self.assertTrue(True)

suite = unittest.TestLoader().loadTestsFromTestCase(TestDbasetoolsFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)

