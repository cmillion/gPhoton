import imagetools as it
import unittest
import numpy as np
import os

# It would be really great to come up with a way to test all of these
#  without assuming that the database is static. If a correction is
#  made to the database at some point, then almost every one of these
#  tests will fail and need to be rewritten.

class TestImagetoolsFunctions(unittest.TestCase):

	def setUp(self):
		self.skyrange = [1,1.5]
		self.skypos   = [20.,120.]

	def test_imgdims(self):
		self.assertEqual(it.imgdims(self.skyrange),[2400.0, 3600.0])
		self.assertEqual(it.imgdims(self.skyrange,width=2001.1,height=2002.8),[2002.0, 2003.0])

	def test_define_wcs(self):
		wcs = it.define_wcs(self.skypos,self.skyrange)
		self.assertAlmostEqual(wcs.wcs.cdelt[0],-0.00041667)
		self.assertAlmostEqual(wcs.wcs.cdelt[1], 0.00041667)
		self.assertEqual(wcs.wcs.ctype[0],'RA---TAN')
		self.assertEqual(wcs.wcs.ctype[1],'DEC--TAN')
		self.assertEqual(wcs.wcs.crpix[0],1200.5)
		self.assertEqual(wcs.wcs.crpix[1],1800.5)
		self.assertEqual(wcs.wcs.crval[0],20.)
		self.assertEqual(wcs.wcs.crval[1],120.)

	def test_movie_tbl(self):
		# FIXME
		self.assertTrue(True)

	def test_fits_header(self):
		# FIXME
		self.assertTrue(True)

	def test_countmap(self):
		# FIXME
		self.assertTrue(True)

	def test_write_jpeg(self):
		# FIXME
		self.assertTrue(True)

	def test_rrhr(self):
		# FIXME
		self.assertTrue(True)

	def test_backgroundmap(self):
		# FIXME
		self.assertTrue(True)

	def test_movie(self):
		# FIXME
		self.assertTrue(True)

	def test_create_images(self):
		# FIXME
		self.assertTrue(True)

	def test_write_images(self):
		# FIXME
		self.assertTrue(True)

suite = unittest.TestLoader().loadTestsFromTestCase(TestImagetoolsFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)

