import galextools as gt
import unittest
import numpy as np
import os

class TestGalextoolsFunctions(unittest.TestCase):

	def setUp(self):
		self.FUV = 'FUV'
		self.NUV = 'NUV'
		self.bands = ['FUV','NUV']
		self.radii = np.arange(35)/1000.

	def test_apcorrect1(self):
		aper = np.array([1.5,2.3,3.8,6.0,9.0,12.8,17.3,30.,60.,90.])/3600.
		# FUV
		dmag = [1.65,0.96,0.36,0.15,0.1,0.09,0.07,0.06,0.03,0.01]
		for i in np.arange(len(aper)):
			self.assertAlmostEqual(gt.apcorrect1(aper[i],'FUV'),dmag[i])
			if aper[i]<aper[-1]:
				self.assertAlmostEqual(gt.apcorrect1((aper[i]+aper[i+1])/2.,'FUV'),(dmag[i]+dmag[i+1])/2.)
		self.assertAlmostEqual(gt.apcorrect1(aper[-1]+10,'FUV'),0.)
		# NUV
		dmag = [2.09,1.33,0.59,0.23,0.13,0.09,0.07,0.04,-0.00,-0.01]
		for i,a in enumerate(aper[:-1]):
			self.assertAlmostEqual(gt.apcorrect1(aper[i],'NUV'),dmag[i])
			if aper[i]<aper[-2]:
				self.assertAlmostEqual(gt.apcorrect1((aper[i]+aper[i+1])/2.,'NUV'),(dmag[i]+dmag[i+1])/2.)
		self.assertAlmostEqual(gt.apcorrect1(aper[-1],'NUV'),0.)
		self.assertAlmostEqual(gt.apcorrect1(aper[-1]+10,'NUV'),0.)

	def test_apcorrect2(self):
		aper = np.array([1.5,2.3,3.8,6.0,9.0,12.8,17.3])/3600.
		# FUV
		band,dmag = 'FUV',[1.65,0.77,0.2,0.1,0.07,0.05,0.04]
		for i in np.arange(len(aper)):
			self.assertAlmostEqual(gt.apcorrect2(aper[i],band),dmag[i])
			if aper[i]<aper[-1]:
				self.assertAlmostEqual(gt.apcorrect2((aper[i]+aper[i+1])/2.,band),(dmag[i]+dmag[i+1])/2.)

		# NUV
		band,dmag = 'NUV',[1.33,0.62,0.21,0.12,0.08,0.06,0.04]
		for i in np.arange(len(aper)):
			self.assertAlmostEqual(gt.apcorrect2(aper[i],band),dmag[i])
			if aper[i]<aper[-1]:
				self.assertAlmostEqual(gt.apcorrect2((aper[i]+aper[i+1])/2.,band),(dmag[i]+dmag[i+1])/2.)

	def test_photometric_repeatability(self):
		cps,expt=10,200
		self.assertAlmostEqual(gt.photometric_repeatability(cps,expt,'FUV'),-2.5*(np.log10(cps)-np.log10(cps+np.sqrt(cps*expt+0.050*cps*expt*0.050*cps*expt)/expt)))
		self.assertAlmostEqual(gt.photometric_repeatability(cps,expt,'NUV'),-2.5*(np.log10(cps)-np.log10(cps+np.sqrt(cps*expt+0.027*cps*expt*0.027*cps*expt)/expt)))

	def test_detbg(self):
		self.assertAlmostEqual(gt.detbg(0.1,'FUV'),0.1 * 1e-4 / ((1.5/(60*60))**2.))
		self.assertAlmostEqual(gt.detbg(0.1,'NUV'),0.1 * 1e-3 / ((1.5/(60*60))**2.))

	def test_counts2mag(self):
		self.assertAlmostEqual(gt.counts2mag(10,'FUV'),-2.5*np.log10(10)+18.82)
		self.assertAlmostEqual(gt.counts2mag(10,'NUV'),-2.5*np.log10(10)+20.08)

	def test_mag2counts(self):
		self.assertAlmostEqual(gt.mag2counts(16,'FUV'),10.**(-(16-18.82)/2.5))
		self.assertAlmostEqual(gt.mag2counts(16,'NUV'),10.**(-(16-20.08)/2.5))

	def test_counts2flux(self):
		self.assertAlmostEqual(gt.counts2flux(10,'FUV'),1.4e-15 * 10)
		self.assertAlmostEqual(gt.counts2flux(10,'NUV'),2.06e-16* 10)

	def test_deg2pix(self):
		self.assertAlmostEqual(gt.deg2pix(1.5),np.ceil(1.5/0.000416666666666667))
		self.assertAlmostEqual(gt.deg2pix(15), np.ceil(15/0.000416666666666667))
		self.assertAlmostEqual(gt.deg2pix(1.5,CDELT2=0.001),np.ceil(1.5/0.001))
		self.assertAlmostEqual(gt.deg2pix(15, CDELT2=0.001),np.ceil(15/0.001))

	def test_compute_flat_scale(self):
		self.assertTrue(True)
		# FUV
		band = 'FUV'
		flat_correct=-0.0154
                flat_t0=840418779.02
                flat_correct_0=1.2420282
                flat_correct_1=-2.8843099e-10
                flat_correct_2=0.000
		for t in [869777733.995,901049156.995]:
			self.assertAlmostEqual(gt.compute_flat_scale(t,band,verbose=0),flat_correct_0+(flat_correct_1*t)+(flat_correct_2*t)*t)
		# NUV
		band = 'NUV'
		flat_correct=-0.0154
		flat_t0=840418779.02
		flat_correct_0=1.9946352
		flat_correct_1=-1.9679445e-09
		flat_correct_2=9.3025231e-19
		for t in [869777733.995,901049156.995]:
			self.assertAlmostEqual(gt.compute_flat_scale(t,band,verbose=0),flat_correct_0+(flat_correct_1*t)+(flat_correct_2*t)*t)

suite = unittest.TestLoader().loadTestsFromTestCase(TestGalextoolsFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)

