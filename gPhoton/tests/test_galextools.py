import gPhoton.galextools as gt
import unittest
import numpy as np
import os

class TestGalextoolsFunctions(unittest.TestCase):

    def setUp(self):
        self.bands = ['FUV','NUV']

    def test_apcorrect1_FUV(self):
        aper = np.array([1.5,2.3,3.8,6.0,9.0,12.8,17.3,30.,60.,90.])/3600.
        dmag = [1.65,0.96,0.36,0.15,0.1,0.09,0.07,0.06,0.03,0.01]
        for i in np.arange(len(aper)):
            self.assertAlmostEqual(gt.apcorrect1(aper[i],'FUV'),dmag[i])
            if aper[i]<aper[-1]:
                self.assertAlmostEqual(gt.apcorrect1((aper[i]+aper[i+1])/2.,'FUV'),(dmag[i]+dmag[i+1])/2.)
        self.assertAlmostEqual(gt.apcorrect1(aper[-1]+10,'FUV'),0.)

    def test_apcorrect1_NUV(self):
        aper = np.array([1.5,2.3,3.8,6.0,9.0,12.8,17.3,30.,60.,90.])/3600.
        dmag = [2.09,1.33,0.59,0.23,0.13,0.09,0.07,0.04,-0.00,-0.01]
        for i,a in enumerate(aper[:-1]):
            self.assertAlmostEqual(gt.apcorrect1(aper[i],'NUV'),dmag[i])
            if aper[i]<aper[-2]:
                self.assertAlmostEqual(gt.apcorrect1((aper[i]+aper[i+1])/2.,'NUV'),(dmag[i]+dmag[i+1])/2.)
        self.assertAlmostEqual(gt.apcorrect1(aper[-1],'NUV'),0.)
        self.assertAlmostEqual(gt.apcorrect1(aper[-1]+10,'NUV'),0.)

    def test_apcorrect2_FUV(self):
        aper = np.array([1.5,2.3,3.8,6.0,9.0,12.8,17.3])/3600.
        band,dmag = 'FUV',[1.65,0.77,0.2,0.1,0.07,0.05,0.04]
        for i in np.arange(len(aper)):
            self.assertAlmostEqual(gt.apcorrect2(aper[i],band),dmag[i])
            if aper[i]<aper[-1]:
                self.assertAlmostEqual(gt.apcorrect2((aper[i]+aper[i+1])/2.,band),(dmag[i]+dmag[i+1])/2.)

    def test_apcorrect2_NUV(self):
        aper = np.array([1.5,2.3,3.8,6.0,9.0,12.8,17.3])/3600.
        band,dmag = 'NUV',[1.33,0.62,0.21,0.12,0.08,0.06,0.04]
        for i in np.arange(len(aper)):
            self.assertAlmostEqual(gt.apcorrect2(aper[i],band),dmag[i])
            if aper[i]<aper[-1]:
                self.assertAlmostEqual(gt.apcorrect2((aper[i]+aper[i+1])/2.,band),(dmag[i]+dmag[i+1])/2.)

    def test_photometric_repeatability_FUV(self):
        cps,expt=10,200
        self.assertAlmostEqual(gt.photometric_repeatability(cps,expt,'FUV'),-2.5*(np.log10(cps)-np.log10(cps+np.sqrt(cps*expt+0.050*cps*expt*0.050*cps*expt)/expt)))

    def test_photometric_repeatability_NUV(self):
        cps,expt=10,200
        self.assertAlmostEqual(gt.photometric_repeatability(cps,expt,'NUV'),-2.5*(np.log10(cps)-np.log10(cps+np.sqrt(cps*expt+0.027*cps*expt*0.027*cps*expt)/expt)))

    def test_detbg_FUV(self):
        self.assertAlmostEqual(gt.detbg(0.1,'FUV'),0.1 * 1e-4 / ((1.5/(60*60))**2.))

    def test_detbg_NUV(self):
        self.assertAlmostEqual(gt.detbg(0.1,'NUV'),0.1 * 1e-3 / ((1.5/(60*60))**2.))

    def test_counts2mag_FUV(self):
        self.assertAlmostEqual(gt.counts2mag(10,'FUV'),-2.5*np.log10(10)+18.82)

    def test_counts2mag_NUV(self):
        self.assertAlmostEqual(gt.counts2mag(10,'NUV'),-2.5*np.log10(10)+20.08)

    def test_mag2counts_FUV(self):
        self.assertAlmostEqual(gt.mag2counts(16,'FUV'),10.**(-(16-18.82)/2.5))

    def test_mag2counts_NUV(self):
        self.assertAlmostEqual(gt.mag2counts(16,'NUV'),10.**(-(16-20.08)/2.5))

    def test_counts2flux_FUV(self):
        self.assertAlmostEqual(gt.counts2flux(10,'FUV'),1.4e-15 * 10)

    def test_counts2flux_NUV(self):
        self.assertAlmostEqual(gt.counts2flux(10,'NUV'),2.06e-16* 10)

    def test_deg2pix_equator(self):
        skypos,skyrange=[0.,0.],[1.,1.]
        pix = gt.deg2pix(skypos,skyrange)
        self.assertAlmostEqual(pix[0],2400.)
        self.assertAlmostEqual(pix[1],2400.)

    def test_deg2pix_pole(self):
        skypos,skyrange=[90.,90.],[1.,1.]
        pix = gt.deg2pix(skypos,skyrange)
        self.assertAlmostEqual(pix[0],2400.)
        self.assertAlmostEqual(pix[1],21.)

    def test_compute_flat_scale_FUV(self):
        band = 'FUV'
        ts=[869777733.995,901049156.995]
        scales=[0.991157347104,0.999816178202]
        for t,scl in zip(ts,scales):
            self.assertAlmostEqual(gt.compute_flat_scale(t,band,verbose=0),scl)

    def test_compute_flat_scale_NUV(self):
        band = 'NUV'
        ts=[869777733.995,901049156.995]
        scales=[0.986709143129,0.994262914909]
        for t,scl in zip(ts,scales):
            self.assertAlmostEqual(gt.compute_flat_scale(t,band,verbose=0),scl)

    def test_compute_flat_scale_FUV_array(self):
        band = 'FUV'
        ts=[869777733.995,901049156.995]
        scales=[0.991157347104,0.999816178202]
        out = gt.compute_flat_scale(ts,band,verbose=0)
        self.assertEqual(len(ts),len(out))
        self.assertAlmostEqual(out[0],scales[0])
        self.assertAlmostEqual(out[0],scales[0])


    def test_compute_flat_scale_NUV_array(self):
        band = 'NUV'
        ts=[869777733.995,901049156.995]
        scales=[0.986709143129,0.994262914909]
        out = gt.compute_flat_scale(ts,band,verbose=0)
        self.assertEqual(len(ts),len(out))
        self.assertAlmostEqual(out[0],scales[0])
        self.assertAlmostEqual(out[0],scales[0])

suite = unittest.TestLoader().loadTestsFromTestCase(TestGalextoolsFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)
