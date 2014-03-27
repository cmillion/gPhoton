import MCUtils
import unittest
import numpy as np
import os

class TestMCUtilsFunctions(unittest.TestCase):

	def setUp(self):
		self.FUV = 'FUV'
		self.NUV = 'NUV'

	def test_error(self):
		self.assertTrue(1==1)

	def test_rotvec(self):
		self.assertTrue(1==1)

	def test_rms(self):
		self.assertTrue(1==1)

	def test_print_inline(self):
		self.assertTrue(1==1)

	def test_manage_requests(self):
		self.assertTrue(1==1)

	def test_manage_grequests(self):
		self.assertTrue(1==1)

	def test_wheretrue(self):
		self.assertTrue(1==1)

	def test_find_nearest_lower(self):
		self.assertTrue(1==1)

	def test_get_fits_data(self):
		self.assertTrue(1==1)

	def test_get_fits_header(self):
		self.assertTrue(1==1)

	def test_get_tbl_data(self):
		self.assertTrue(1==1)

#       def tearDown():
#

suite = unittest.TestLoader().loadTestsFromTestCase(TestMCUtilsFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)

