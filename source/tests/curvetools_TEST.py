import curvetools as ct
import unittest
import numpy as np
import os

# It would be really great to come up with a way to test all of these
#  without assuming that the database is static. If a correction is
#  made to the database at some point, then almost every one of these
#  tests will fail and need to be rewritten.

class TestCurvetoolsFunctions(unittest.TestCase):

	def setUp(self):
		self.foo = 1



suite = unittest.TestLoader().loadTestsFromTestCase(TestCurvetoolsFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)

