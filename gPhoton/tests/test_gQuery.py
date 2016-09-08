from __future__ import absolute_import, division, print_function
import unittest
import gPhoton.gQuery as gq
import gPhoton.CalUtils as cu
from gPhoton import time_id

class TestGQueryFunctions(unittest.TestCase):
    def setUp(self):
        self.bands = ['NUV','FUV']
        self.NUV = 'NUV'
        self.FUV = 'FUV'
        self.ra0 = 176.919525856
        self.dec0 = 0.255696872807
        self.t0 = 766525332.995
        self.t1 = 866526576.995
        self.radius = 0.004
        self.eclipse = 23456
        self.maglimit = 30.
        self.detsize = 1.25
        self.xr = [200,400]
        self.yr = [300,500]
        self.tscale = 1000.
        self.aspum=68.754932/1000.
        self.baseURL = 'https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query='
        self.baseDB = 'GPFCore.dbo'
        self.MCATDB = 'GR6Plus7.dbo'
        self.formatURL = ' -- '+str(time_id)+'&format=extjs'

    def test_baseURL(self):
        self.assertEqual(gq.baseURL,self.baseURL)

    def test_formatURL(self):
        self.assertEqual(gq.formatURL,self.formatURL)

    def test_mcat_sources(self):
        band = self.NUV
        bandflag = 1 if band=='NUV' else 2
        query = (str(self.baseURL)+
            'select ra, dec, nuv_mag, fuv_mag, fov_radius, nuv_skybg, fuv_skybg,'
            ' nuv_fwhm_world, fuv_fwhm_world, fuv_mag_aper_1, fuv_mag_aper_2,'
            ' fuv_mag_aper_3, fuv_mag_aper_4, fuv_mag_aper_5, fuv_mag_aper_6,'
            ' fuv_mag_aper_7, nuv_mag_aper_1, nuv_mag_aper_2, nuv_mag_aper_3,'
            ' nuv_mag_aper_4, nuv_mag_aper_5, nuv_mag_aper_6, nuv_mag_aper_7'
            ' fuv_magerr_aper_1, fuv_magerr_aper_2, fuv_magerr_aper_3,'
            ' fuv_magerr_aper_4, fuv_magerr_aper_5, fuv_magerr_aper_6,'
            ' fuv_magerr_aper_7, nuv_magerr_aper_1, nuv_magerr_aper_2,'
            ' nuv_magerr_aper_3, nuv_magerr_aper_4, nuv_magerr_aper_5,'
            ' nuv_magerr_aper_6, nuv_magerr_aper_7'
            ' from '+str(self.MCATDB)+'.photoobjall as p inner join '+str(self.MCATDB)+
            '.photoextract as pe on p.photoextractid=pe.photoextractid inner join '+
            str(self.MCATDB)+'.fgetnearbyobjeq('+repr(float(self.ra0))+', '+
            repr(float(self.dec0))+', '+
            str(self.radius*60.)+') as nb on p.objid=nb.objid and (band=3 or band='+
            str(bandflag)+') and '+str(band)+'_mag<'+str(self.maglimit)+
            str(self.formatURL))
        self.assertEqual(gq.mcat_sources(self.NUV,self.ra0,self.dec0,self.radius,maglimit=self.maglimit),query)

    def test_exposure_ranges(self):
        query = "https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select distinct time from GPFCore.dbo.fGetNearbyAspectEq(176.919525856,0.255696872807,((1.25/2.0)*60.0),766525332995,866526576996) where band='NUV' or band='FUV/NUV' order by time"+self.formatURL
        self.assertEqual(gq.exposure_ranges(self.NUV,self.ra0,self.dec0,t0=self.t0,t1=self.t1,detsize=self.detsize),query)

    def test_exposure_range(self):
        self.assertEqual(gq.exposure_range(self.NUV,self.ra0,self.dec0,t0=self.t0,t1=self.t1),"https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select startTimeRange, endTimeRange from GPFCore.dbo.fGetTimeRanges(766525332,866526576,176.919525856,0.255696872807) where band='NUV'"+self.formatURL)

    def test_aperture(self):
        self.assertEqual(gq.aperture(self.NUV,self.ra0,self.dec0,self.t0,self.t1,self.radius),'https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=Select  sum(photonCount) from GPFCore.dbo.fGetNearbyObjEqCountNUV(176.919525856,0.255696872807,0.004,766525332995,866526576995,0)'+self.formatURL)

    def test_deadtime1(self):
        self.assertEqual(gq.deadtime1(self.NUV,self.t0,self.t1),'https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select count(*) from GPFCore.dbo.NUVPhotonsV where time >= 766525332995 and time < 866526576995'+self.formatURL)

    def test_deadtime2(self):
        self.assertEqual(gq.deadtime2(self.NUV,self.t0,self.t1),'https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select count(*) from GPFCore.dbo.NUVPhotonsNULLV where time >= 766525332995 and time < 866526576995'+self.formatURL)

    def test_deadtime(self):
        self.assertEqual(gq.deadtime(self.NUV,self.t0,self.t1),'https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select sum(dt)*0.000005714285714285714511150895 / (866526576.995-766525332.995) from(select count(*) as dt from GPFCore.dbo.NUVPhotonsNULLV where time >= 766525332995 and time < 866526576995 union all select count(*) as dt from GPFCore.dbo.NUVPhotonsV where time >= 766525332995 and time < 866526576995) x'+self.formatURL)

    def test_boxcount(self):
        self.assertEqual(gq.boxcount(self.NUV,self.t0,self.t1,self.xr,self.yr),'https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select count(*) from GPFCore.dbo.NUVPhotonsNULLV where time >= 766525332995 and time < 866526576995 and x >= 200 and x < 400 and y >= 300 and y < 500'+self.formatURL)

    def test_stimcount(self):
        margin=[90.01,90.01]
        avgstim = cu.avg_stimpos(self.NUV,self.eclipse)
        query = ('{baseURL}select count(*) from {baseDB}.{band}PhotonsNULLV '+
            'where time >= {t0} and time < {t1} and ('+
            '((x >= {x10} and x < {x11}) and (y >= {y10} and y < {y11})) or '+
            '((x >= {x20} and x < {x21}) and (y >= {y20} and y < {y21})) or '+
            '((x >= {x30} and x < {x31}) and (y >= {y30} and y < {y31})) or '+
            '((x >= {x40} and x < {x41}) and (y >= {y40} and y < {y41}))'+
            '){formatURL}').format(baseURL=self.baseURL,
                baseDB=self.baseDB, band=self.NUV,
                t0=str(long(self.t0*self.tscale)),
                t1=str(long(self.t1*self.tscale)),
                x10=(avgstim['x1']-margin[0])/self.aspum,
                x11=(avgstim['x1']+margin[0])/self.aspum,
                y10=(avgstim['y1']-margin[1])/self.aspum,
                y11=(avgstim['y1']+margin[1])/self.aspum,
                x20=(avgstim['x2']-margin[0])/self.aspum,
                x21=(avgstim['x2']+margin[0])/self.aspum,
                y20=(avgstim['y2']-margin[1])/self.aspum,
                y21=(avgstim['y2']+margin[1])/self.aspum,
                x30=(avgstim['x3']-margin[0])/self.aspum,
                x31=(avgstim['x3']+margin[0])/self.aspum,
                y30=(avgstim['y3']-margin[1])/self.aspum,
                y31=(avgstim['y3']+margin[1])/self.aspum,
                x40=(avgstim['x4']-margin[0])/self.aspum,
                x41=(avgstim['x4']+margin[0])/self.aspum,
                y40=(avgstim['y4']-margin[1])/self.aspum,
                y41=(avgstim['y4']+margin[1])/self.aspum,
                formatURL=self.formatURL)
        self.assertEqual(gq.stimcount(self.NUV,self.t0,self.t1,eclipse=self.eclipse),query)

    def test_boxcentroid(self):
        self.assertEqual(gq.boxcentroid(self.NUV,self.t0,self.t1,self.xr,self.yr),'https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select avg(x), avg(y) from GPFCore.dbo.NUVPhotonsNULLV where time >= 766525332995 and time < 866526576995 and x >= 200 and x < 400 and y >= 300 and y < 500'+self.formatURL)

    def test_boxtimes(self):
        self.assertEqual(gq.boxtimes(self.NUV,self.t0,self.t1,self.xr,self.yr),'https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select time from GPFCore.dbo.NUVPhotonsNULLV where time >= 766525332995 and time < 866526576995 and x >= 200 and x < 400 and y >= 300 and y < 500'+self.formatURL)

    def test_allphotons(self):
        self.assertEqual(gq.allphotons(self.NUV,self.ra0,self.dec0,self.t0,self.t1,self.radius),'https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select time,ra,dec,xi,eta,x,y,flag from GPFCore.dbo.fGetNearbyObjEqNUVAllColumns(176.919525856,0.255696872807,0.004,766525332995,866526576995,0)'+self.formatURL)

    def test_shutter(self):
        self.assertEqual(gq.shutter(self.NUV,self.t0,self.t1),'https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select shutter*0.05 from GPFCore.dbo.fGetNUVShutter(766525332995,866526576995)'+self.formatURL)

    def test_aspect(self):
        self.assertEqual(gq.aspect(self.t0,self.t1),'https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select eclipse, filename, time, ra, dec, twist, flag, ra0, dec0, twist0 from aspect where time >= 766525332995 and time < 866526576995 order by time'+self.formatURL)

    def test_aspect_ecl(self):
        self.assertEqual(gq.aspect_ecl(self.eclipse),'https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select eclipse, filename, time, ra, dec, twist, flag, ra0, dec0, twist0 from aspect where eclipse=23456 order by time'+self.formatURL)

    def test_box(self):
        self.assertEqual(gq.box(self.NUV,self.ra0,self.dec0,self.t0,self.t1,self.radius),'https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select time,ra,dec from GPFCore.dbo.NUVPhotonsV where time >= 766525332995 and time < 866526576995 and ra >= 176.91552585600002 and ra < 176.923525856 and dec >= 0.251696872807 and dec < 0.259696872807 and flag=0'+self.formatURL)

suite = unittest.TestLoader().loadTestsFromTestCase(TestGQueryFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)
