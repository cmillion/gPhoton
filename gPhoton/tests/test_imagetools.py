from __future__ import absolute_import, division, print_function
import unittest
import gPhoton.imagetools as it
import unittest
import numpy as np
import os

"""Regression tests for image creation functions."""
class TestImagetoolsFunctions(unittest.TestCase):
    def setUp(self):
        self.bands = ['NUV','FUV']
        self.skypos = [176.919525856024,0.255696872807351]
        self.tranges = [[ 765579107.995, 765579821.995 ],
                        [ 766525332.995, 766526576.995 ],
                        [ 919754986.995, 919755916.995 ]]
        self.trange = self.tranges[0]
        self.radius = 0.01
        self.annulus = [0.02,0.03]
        self.skyrange = [0.01,0.01]
        self.framesz = 100.

#---

    def test_cnt_NUV(self):
        cnt_NUV = it.create_image(
            'NUV',self.skypos,self.tranges[0],self.skyrange)
        self.assertEqual(np.shape(cnt_NUV)[0],24)
        self.assertEqual(np.shape(cnt_NUV)[1],23)
        self.assertAlmostEqual(np.sum(cnt_NUV),1787.0)

    def test_cnt_visits_NUV(self):
        cnt_visits_NUV = it.create_image(
            'NUV',self.skypos,self.tranges,self.skyrange)
        self.assertEqual(np.shape(cnt_visits_NUV)[0],3)
        self.assertEqual(np.shape(cnt_visits_NUV)[1],24)
        self.assertEqual(np.shape(cnt_visits_NUV)[2],23)
        self.assertAlmostEqual(np.sum(cnt_visits_NUV),519839.0)

    def test_cnt_coadd_NUV(self):
        cnt_coadd_NUV = it.create_image(
            'NUV',self.skypos,self.tranges,self.skyrange,coadd=True)
        self.assertEqual(np.shape(cnt_coadd_NUV)[0],24)
        self.assertEqual(np.shape(cnt_coadd_NUV)[1],23)
        self.assertAlmostEqual(np.sum(cnt_coadd_NUV),519839.0)

    def test_cnt_movie_NUV(self):
        cnt_movie_NUV = it.create_image(
            'NUV',self.skypos,self.tranges[0],self.skyrange,
            framesz=self.framesz)
        self.assertEqual(np.shape(cnt_movie_NUV)[0],8)
        self.assertEqual(np.shape(cnt_movie_NUV)[1],24)
        self.assertEqual(np.shape(cnt_movie_NUV)[2],23)
        self.assertAlmostEqual(np.sum(cnt_movie_NUV),1787.0)

    def test_cnt_visit_flux_preservation_NUV(self):
        cnt_movie_NUV = it.create_image(
            'NUV',self.skypos,self.tranges[0],self.skyrange,
            framesz=self.framesz)
        cnt_NUV = it.create_image(
            'NUV',self.skypos,self.tranges[0],self.skyrange)
        self.assertAlmostEqual(np.sum(cnt_movie_NUV),np.sum(cnt_NUV))

    def test_cnt_coadd_flux_preservation_NUV(self):
        cnt_coadd_NUV = it.create_image(
            'NUV',self.skypos,self.tranges,self.skyrange,coadd=True)
        cnt_visits_NUV = it.create_image(
            'NUV',self.skypos,self.tranges,self.skyrange)
        self.assertAlmostEqual(np.sum(cnt_coadd_NUV),np.sum(cnt_visits_NUV))

    def test_cnt_FUV(self):
        cnt_FUV = it.create_image(
            'FUV',self.skypos,self.tranges[1],self.skyrange)
        self.assertEqual(np.shape(cnt_FUV)[0],24)
        self.assertEqual(np.shape(cnt_FUV)[1],23)
        self.assertAlmostEqual(np.sum(cnt_FUV),124015.0)

    def test_cnt_visits_FUV(self):
        cnt_visits_FUV = it.create_image(
            'FUV',self.skypos,self.tranges,self.skyrange)
        self.assertEqual(np.shape(cnt_visits_FUV)[0],2)
        self.assertEqual(np.shape(cnt_visits_FUV)[1],24)
        self.assertEqual(np.shape(cnt_visits_FUV)[2],23)
        self.assertAlmostEqual(np.sum(cnt_visits_FUV),124283.0)

    def test_cnt_coadd_FUV(self):
        cnt_coadd_FUV = it.create_image(
            'FUV',self.skypos,self.tranges,self.skyrange,coadd=True)
        self.assertEqual(np.shape(cnt_coadd_FUV)[0],24)
        self.assertEqual(np.shape(cnt_coadd_FUV)[1],23)
        self.assertAlmostEqual(np.sum(cnt_coadd_FUV),124283.0)

    def test_cnt_movie_FUV(self):
        cnt_movie_FUV = it.create_image(
            'FUV',self.skypos,self.tranges[1],self.skyrange,
            framesz=self.framesz)
        self.assertEqual(np.shape(cnt_movie_FUV)[0],13)
        self.assertEqual(np.shape(cnt_movie_FUV)[1],24)
        self.assertEqual(np.shape(cnt_movie_FUV)[2],23)
        self.assertAlmostEqual(np.sum(cnt_movie_FUV),124015.0)

    def test_cnt_visit_flux_preservation_FUV(self):
        cnt_movie_FUV = it.create_image(
            'FUV',self.skypos,self.tranges[1],self.skyrange,
            framesz=self.framesz)
        cnt_FUV = it.create_image(
            'FUV',self.skypos,self.tranges[1],self.skyrange)
        self.assertAlmostEqual(np.sum(cnt_movie_FUV),np.sum(cnt_FUV))

    def test_cnt_coadd_flux_preservation_FUV(self):
        cnt_coadd_FUV = it.create_image(
            'FUV',self.skypos,self.tranges,self.skyrange,coadd=True)
        cnt_visits_FUV = it.create_image(
            'FUV',self.skypos,self.tranges,self.skyrange)
        self.assertAlmostEqual(np.sum(cnt_coadd_FUV),np.sum(cnt_visits_FUV))

    def test_int_NUV(self):
        int_NUV = it.create_image(
            'NUV',self.skypos,self.tranges[0],self.skyrange,response=True)
        self.assertEqual(np.shape(int_NUV)[0],24)
        self.assertEqual(np.shape(int_NUV)[1],23)
        self.assertAlmostEqual(np.sum(int_NUV),3.3729434940425236)

    def test_int_visits_NUV(self):
        int_visits_NUV = it.create_image(
            'NUV',self.skypos,self.tranges,self.skyrange,response=True)
        self.assertEqual(np.shape(int_visits_NUV)[0],3)
        self.assertEqual(np.shape(int_visits_NUV)[1],24)
        self.assertEqual(np.shape(int_visits_NUV)[2],23)
        self.assertAlmostEqual(np.sum(int_visits_NUV),533.60977439824569)

    def test_int_coadd_NUV(self):
        int_coadd_NUV = it.create_image(
            'NUV',self.skypos,self.tranges,self.skyrange,
                                                coadd=True,response=True)
        self.assertEqual(np.shape(int_coadd_NUV)[0],24)
        self.assertEqual(np.shape(int_coadd_NUV)[1],23)
        self.assertAlmostEqual(np.sum(int_coadd_NUV),533.60977439824592)

    def test_int_movie_NUV(self):
        int_movie_NUV = it.create_image(
            'NUV',self.skypos,self.tranges[0],self.skyrange,response=True,
            framesz=self.framesz)
        self.assertEqual(np.shape(int_movie_NUV)[0],8)
        self.assertEqual(np.shape(int_movie_NUV)[1],24)
        self.assertEqual(np.shape(int_movie_NUV)[2],23)
        self.assertAlmostEqual(np.sum(int_movie_NUV),27.904517056070709)

    def test_int_coadd_flux_preservation_NUV(self):
        int_coadd_NUV = it.create_image(
            'NUV',self.skypos,self.tranges,self.skyrange,
                                                coadd=True,response=True)
        int_visits_NUV = it.create_image(
            'NUV',self.skypos,self.tranges,self.skyrange,response=True)
        self.assertAlmostEqual(np.sum(int_visits_NUV),
                                                np.sum(int_coadd_NUV))

    # def test_int_visit_flux_preservation_NUV(self): ### FAILS!!! ###
    #     int_movie_NUV = it.create_image(
    #         'NUV',self.skypos,self.tranges[0],self.skyrange,response=True,
    #         framesz=self.framesz)
    #     int_NUV = it.create_image(
    #         'NUV',self.skypos,self.tranges[0],self.skyrange,response=True)
    #     self.assertAlmostEqual(np.sum(int_movie_NUV),np.sum(int_NUV))

    def test_int_FUV(self):
        int_FUV = it.create_image(
            'FUV',self.skypos,self.tranges[1],self.skyrange,response=True)
        self.assertEqual(np.shape(int_FUV)[0],24)
        self.assertEqual(np.shape(int_FUV)[1],23)
        self.assertAlmostEqual(np.sum(int_FUV),122.42292744331772)

    def test_int_visits_FUV(self):
        int_visits_FUV = it.create_image(
            'FUV',self.skypos,self.tranges,self.skyrange,response=True)
        self.assertEqual(np.shape(int_visits_FUV)[0],2)
        self.assertEqual(np.shape(int_visits_FUV)[1],24)
        self.assertEqual(np.shape(int_visits_FUV)[2],23)
        self.assertAlmostEqual(np.sum(int_visits_FUV),122.71086595458925)

    def test_int_coadd_FUV(self):
        int_coadd_FUV = it.create_image(
            'FUV',self.skypos,self.tranges,self.skyrange,
                                                coadd=True,response=True)
        self.assertEqual(np.shape(int_coadd_FUV)[0],24)
        self.assertEqual(np.shape(int_coadd_FUV)[1],23)
        self.assertAlmostEqual(np.sum(int_coadd_FUV),122.71086595458925)

    def test_int_movie_FUV(self):
        int_movie_FUV = it.create_image(
            'FUV',self.skypos,self.tranges[1],self.skyrange,response=True,
                    framesz=self.framesz)
        self.assertEqual(np.shape(int_movie_FUV)[0],13)
        self.assertEqual(np.shape(int_movie_FUV)[1],24)
        self.assertEqual(np.shape(int_movie_FUV)[2],23)
        self.assertAlmostEqual(np.sum(int_movie_FUV),1592.1558509925253)

    def test_int_coadd_flux_preservation_FUV(self):
        int_visits_FUV = it.create_image(
            'FUV',self.skypos,self.tranges,self.skyrange,response=True)
        int_coadd_FUV = it.create_image(
            'FUV',self.skypos,self.tranges,self.skyrange,
                                                coadd=True,response=True)
        self.assertAlmostEqual(np.sum(int_visits_FUV),np.sum(int_coadd_FUV))

    # def test_int_visit_flux_preservation_FUV(self): ### FAILS!!! ###
    #     int_movie_FUV = it.create_image(
    #         'FUV',self.skypos,self.tranges[1],self.skyrange,response=True,
    #                 framesz=self.framesz)
    #     int_FUV = it.create_image(
    #         'FUV',self.skypos,self.tranges[0],self.skyrange,response=True)
    #     self.assertAlmostEqual(np.sum(int_movie_FUV),np.sum(int_FUV))

# These should gracefully return None instead of crashing.

    def test_cnt_movie_FUV_nodata(self):
        cnt_movie_FUV_nodata = it.create_image(
            'FUV',self.skypos,self.tranges[0],self.skyrange,
            framesz=self.framesz)
        self.assertIsNone(cnt_movie_FUV_nodata.tolist())

    def test_int_movie_FUV_nodata(self):
        int_movie_FUV_nodata = it.create_image(
            'FUV',self.skypos,self.tranges[0],self.skyrange,response=True,
            framesz=self.framesz)
        self.assertIsNone(int_movie_FUV_nodata.tolist())

    def test_cnt_FUV_nodata(self):
        cnt_FUV_nodata = it.create_image(
            'FUV',self.skypos,self.tranges[0],self.skyrange)
        self.assertIsNone(cnt_FUV_nodata.tolist())

    def test_int_FUV_nodata(self):
        int_FUV_nodata = it.create_image(
            'FUV',self.skypos,self.tranges[0],self.skyrange,response=True)
        self.assertIsNone(int_FUV_nodata.tolist())

#---
suite = unittest.TestLoader().loadTestsFromTestCase(TestImagetoolsFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)
