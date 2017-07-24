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
        for i, frame in enumerate(int_visits_NUV):
            self.assertAlmostEqual(
                [3.37294349404, 527.037522816, 3.19930808804][i],
                np.sum(frame))

    def test_int_coadd_NUV(self):
        int_coadd_NUV = it.create_image(
            'NUV',self.skypos,self.tranges,self.skyrange,
                                                coadd=True,response=True)
        self.assertEqual(np.shape(int_coadd_NUV)[0],24)
        self.assertEqual(np.shape(int_coadd_NUV)[1],23)
        self.assertAlmostEqual(np.sum(int_coadd_NUV),228.12301821935804)

    def test_int_movie_NUV(self):
        int_movie_NUV = it.create_image(
            'NUV',self.skypos,self.tranges[0],self.skyrange,response=True,
            framesz=self.framesz)
        self.assertEqual(np.shape(int_movie_NUV)[0],8)
        self.assertEqual(np.shape(int_movie_NUV)[1],24)
        self.assertEqual(np.shape(int_movie_NUV)[2],23)
        for i,frame in enumerate(int_movie_NUV):
            self.assertAlmostEqual(
                [3.5552140888, 3.27573841449, 3.46923681959, 3.35124651346,
                 3.06126920207, 3.24649727587, 3.50245417683, 4.44126218936][i],
                np.sum(frame))

    def test_int_coadd_flux_preservation_NUV(self):
        """ Intensity maps are basically countrate maps, which is to say that
        they are averaged in time. The sum of the averages is not the average
        of the sums, so these should not be equal unless exposure time
        correction is not happening and we're accidentally creating a "flux"
        map.
        """
        int_coadd_NUV = it.create_image(
            'NUV',self.skypos,self.tranges,self.skyrange,
                                                coadd=True,response=True)
        int_visits_NUV = it.create_image(
            'NUV',self.skypos,self.tranges,self.skyrange,response=True)
        self.assertNotEqual(np.sum(int_visits_NUV),
                                                np.sum(int_coadd_NUV))

    def test_int_visit_flux_preservation_NUV(self):
        """ Intensity maps are basically countrate maps, which is to say that
        they are averaged in time. The sum of the averages is not the average
        of the sums, so these should not be equal unless exposure time
        correction is not happening and we're accidentally creating a "flux"
        map.
        """
        int_movie_NUV = it.create_image(
            'NUV',self.skypos,self.tranges[0],self.skyrange,response=True,
            framesz=self.framesz)
        int_NUV = it.create_image(
            'NUV',self.skypos,self.tranges[0],self.skyrange,response=True)
        self.assertNotEqual(np.sum(int_movie_NUV),np.sum(int_NUV))

    def test_int_FUV(self):
        int_FUV = it.create_image(
            'FUV',self.skypos,self.tranges[1],self.skyrange,response=True)
        self.assertEqual(np.shape(int_FUV)[0],24)
        self.assertEqual(np.shape(int_FUV)[1],23)
        self.assertAlmostEqual(np.sum(int_FUV),122.42292744331783)

    def test_int_visits_FUV(self):
        int_visits_FUV = it.create_image(
            'FUV',self.skypos,self.tranges,self.skyrange,response=True)
        self.assertEqual(np.shape(int_visits_FUV)[0],2)
        self.assertEqual(np.shape(int_visits_FUV)[1],24)
        self.assertEqual(np.shape(int_visits_FUV)[2],23)
        for i,frame in enumerate(int_visits_FUV):
            self.assertAlmostEqual([122.422927443, 0.287938511271][i],
                np.sum(frame))

    def test_int_coadd_FUV(self):
        int_coadd_FUV = it.create_image(
            'FUV',self.skypos,self.tranges,self.skyrange,
                                                coadd=True,response=True)
        self.assertEqual(np.shape(int_coadd_FUV)[0],24)
        self.assertEqual(np.shape(int_coadd_FUV)[1],23)
        self.assertAlmostEqual(np.sum(int_coadd_FUV),70.043094251503163)

    def test_int_movie_FUV(self):
        int_movie_FUV = it.create_image(
            'FUV',self.skypos,self.tranges[1],self.skyrange,response=True,
                    framesz=self.framesz)
        self.assertEqual(np.shape(int_movie_FUV)[0],13)
        self.assertEqual(np.shape(int_movie_FUV)[1],24)
        self.assertEqual(np.shape(int_movie_FUV)[2],23)
        for i,frame in enumerate(int_movie_FUV):
            self.assertAlmostEqual(
                [0.497924179493, 0.847007693848, 0.80058446059, 0.972137734434,
                 1.3665057898, 5.39927107707, 171.572488917, 195.807659308,
                 190.440959947, 455.239978231, 282.896126611, 164.741134199,
                 121.560357866][i], np.sum(frame))

    def test_int_coadd_flux_preservation_FUV(self):
        """ Intensity maps are basically countrate maps, which is to say that
        they are averaged in time. The sum of the averages is not the average
        of the sums, so these should not be equal unless exposure time
        correction is not happening and we're accidentally creating a "flux"
        map.
        """
        int_visits_FUV = it.create_image(
            'FUV',self.skypos,self.tranges,self.skyrange,response=True)
        int_coadd_FUV = it.create_image(
            'FUV',self.skypos,self.tranges,self.skyrange,
                                                coadd=True,response=True)
        self.assertNotEqual(np.sum(int_visits_FUV),np.sum(int_coadd_FUV))

    def test_int_visit_flux_preservation_FUV(self): ### FAILS!!! ###
        """ Intensity maps are basically countrate maps, which is to say that
        they are averaged in time. The sum of the averages is not the average
        of the sums, so these should not be equal unless exposure time
        correction is not happening and we're accidentally creating a "flux"
        map.
        """
        int_movie_FUV = it.create_image(
            'FUV',self.skypos,self.tranges[1],self.skyrange,response=True,
                    framesz=self.framesz)
        int_FUV = it.create_image(
            'FUV',self.skypos,self.tranges[1],self.skyrange,response=True)
        self.assertNotEqual(np.sum(int_movie_FUV),np.sum(int_FUV))

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
