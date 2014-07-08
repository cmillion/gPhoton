import FileUtils
import unittest
import numpy as np
import os

class TestFileUtilsFunctions(unittest.TestCase):

	def setUp(self):
		self.NUV = 'NUV'
		self.FUV = 'FUV'
		# a random test eclipse
		self.eclipse = 31000
		# an eclipse that is before the CSP
		self.preCSP  = 31000
		# an eclipse that is after the CSP
		self.postCSP = 40000
		self.detsize = 1.5
		self.imsz    = 3840.
		self.pixsz   = 0.00016666667
		self.npixx   = 3840.
		self.npixy   = 3840.
		self.fill    = self.detsize/(self.npixx*self.pixsz)
		self.raw6_fuv= '../e31000/MISWZN11_12494_0315_0002-fd-raw6.fits'
		self.raw6_nuv= '../e31000/MISWZN11_12494_0315_0002-nd-raw6.fits'
		self.aspfile = '../e31000/MISWZN11_12494_0315_0002-asprta.fits'
		self.ssd  = 'testSSD.tbl'
		self.calpath = '../cal/'
		self.pltscl = 68.754932
		self.aspum  = self.pltscl/1000.

	def test_find_band(self):
		self.assertEqual(FileUtils.find_band(self.raw6_fuv),self.FUV)
		self.assertEqual(FileUtils.find_band(self.raw6_nuv),self.NUV)

	def test_load_aspect(self):
		# FIXME
		self.assertTrue(True)

	def test_web_query_aspect(self):
		asp = FileUtils.web_query_aspect(self.eclipse)
		self.assertEqual(len(asp),6)
		self.assertEqual(len(asp[0]),1445)
		self.assertEqual(asp[5][100],9.0)
		self.assertEqual(asp[3][500],918992208.995)

	def test_wiggle_filename(self):
		self.assertEqual(FileUtils.wiggle_filenames(self.FUV,self.calpath),{'y': '../cal/FUV_wiggle_y.fits', 'x': '../cal/FUV_wiggle_x.fits'})
		self.assertEqual(FileUtils.wiggle_filenames(self.NUV,self.calpath),{'y': '../cal/NUV_wiggle_y.fits', 'x': '../cal/NUV_wiggle_x.fits'})

	def test_avgwalk_filenames(self):
		self.assertEqual(FileUtils.avgwalk_filenames(self.FUV,self.calpath),{'y': '../cal/FUV_avgwalk_y.fits', 'x': '../cal/FUV_avgwalk_x.fits'})
		self.assertEqual(FileUtils.avgwalk_filenames(self.NUV,self.calpath),{'y': '../cal/NUV_avgwalk_y.fits', 'x': '../cal/NUV_avgwalk_x.fits'})

	def test_walk_filenames(self):
		self.assertEqual(FileUtils.walk_filenames(self.FUV,self.calpath),{'y': '../cal/FUV_walk_y.fits', 'x': '../cal/FUV_walk_x.fits'})
		self.assertEqual(FileUtils.walk_filenames(self.NUV,self.calpath),{'y': '../cal/NUV_walk_y.fits', 'x': '../cal/NUV_walk_x.fits'})

	def test_linearity_filenames(self):
		self.assertEqual(FileUtils.linearity_filenames(self.FUV,self.calpath),{'y': '../cal/FUV_NLC_y_det2sky.fits', 'x': '../cal/FUV_NLC_x_det2sky.fits'})
		self.assertEqual(FileUtils.linearity_filenames(self.NUV,self.calpath),{'y': '../cal/NUV_NLC_y_det2sky.fits', 'x': '../cal/NUV_NLC_x_det2sky.fits'})

	def test_flat_filenames(self):
		self.assertEqual(FileUtils.flat_filename(self.FUV,self.calpath),'../cal/FUV_flat.fits')
		self.assertEqual(FileUtils.flat_filename(self.NUV,self.calpath),'../cal/NUV_flat.fits')

	def test_distortion_filenames(self):
		# The FUV filename is always the same.
		self.assertEqual(FileUtils.distortion_filenames(self.FUV,self.calpath,self.preCSP,5000.0),{'y': '../cal/fuv_distortion_cube_dy.fits', 'x': '../cal/fuv_distortion_cube_dx.fits'})
		self.assertEqual(FileUtils.distortion_filenames(self.FUV,self.calpath,self.preCSP,5136.5),{'y': '../cal/fuv_distortion_cube_dy.fits', 'x': '../cal/fuv_distortion_cube_dx.fits'})
		self.assertEqual(FileUtils.distortion_filenames(self.FUV,self.calpath,self.preCSP,6000.0),{'y': '../cal/fuv_distortion_cube_dy.fits', 'x': '../cal/fuv_distortion_cube_dx.fits'})
		self.assertEqual(FileUtils.distortion_filenames(self.FUV,self.calpath,self.postCSP,5000.0),{'y': '../cal/fuv_distortion_cube_dy.fits', 'x': '../cal/fuv_distortion_cube_dx.fits'})
		self.assertEqual(FileUtils.distortion_filenames(self.FUV,self.calpath,self.postCSP,5136.5),{'y': '../cal/fuv_distortion_cube_dy.fits', 'x': '../cal/fuv_distortion_cube_dx.fits'})
		self.assertEqual(FileUtils.distortion_filenames(self.FUV,self.calpath,self.postCSP,6000.0),{'y': '../cal/fuv_distortion_cube_dy.fits', 'x': '../cal/fuv_distortion_cube_dx.fits'})
		# The NUV filename is the same before the CSP
		self.assertEqual(FileUtils.distortion_filenames(self.NUV,self.calpath,self.preCSP,5000.0),{'y': '../cal/nuv_distortion_cube_dy.fits', 'x': '../cal/nuv_distortion_cube_dx.fits'})
		self.assertEqual(FileUtils.distortion_filenames(self.NUV,self.calpath,self.preCSP,5136.5),{'y': '../cal/nuv_distortion_cube_dy.fits', 'x': '../cal/nuv_distortion_cube_dx.fits'})
		self.assertEqual(FileUtils.distortion_filenames(self.NUV,self.calpath,self.preCSP,6000.0),{'y': '../cal/nuv_distortion_cube_dy.fits', 'x': '../cal/nuv_distortion_cube_dx.fits'})
		# The NUV filename varies by stim separation after the CSP
		self.assertEqual(FileUtils.distortion_filenames(self.NUV,self.calpath,self.postCSP,5000.0),{'y': '../cal/nuv_distortion_cube_dya.fits', 'x': '../cal/nuv_distortion_cube_dxa.fits'})
		self.assertEqual(FileUtils.distortion_filenames(self.NUV,self.calpath,self.postCSP,5136.5),{'y': '../cal/nuv_distortion_cube_dyb.fits', 'x': '../cal/nuv_distortion_cube_dxb.fits'})
		self.assertEqual(FileUtils.distortion_filenames(self.NUV,self.calpath,self.postCSP,6000.0),{'y': '../cal/nuv_distortion_cube_dyc.fits', 'x': '../cal/nuv_distortion_cube_dxc.fits'})

	def test_offset_filenames(self):
		self.assertEqual(FileUtils.offset_filenames(self.calpath),{'y': '../cal/fuv_dy_fdttdc_coef_0.tbl', 'x': '../cal/fuv_dx_fdttdc_coef_0.tbl'})

	def test_mask_filename(self):
		self.assertEqual(FileUtils.mask_filename(self.FUV,self.calpath),'../cal/FUV_mask.fits')
		self.assertEqual(FileUtils.mask_filename(self.NUV,self.calpath),'../cal/NUV_mask.fits')

	def test_create_SSD_filename(self):
		self.assertEqual(FileUtils.create_SSD_filename(self.FUV,self.eclipse),'SSD_fuv_31000.tbl')
		self.assertEqual(FileUtils.create_SSD_filename(self.NUV,self.eclipse),'SSD_nuv_31000.tbl')

suite = unittest.TestLoader().loadTestsFromTestCase(TestFileUtilsFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)


