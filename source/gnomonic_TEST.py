import gnomonic
import unittest
import numpy as np
import os

class TestGnomonicFunctions(unittest.TestCase):

	def setUp(self):
		self.foo = 1

	def test_gnomrev_simple(self):
		self.assertTrue(True)

	def test_gnomfwd_simple(self):
		self.assertTrue(True)


suite = unittest.TestLoader().loadTestsFromTestCase(TestGnomonicFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)

