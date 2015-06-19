import unittest
import gQuery as gq

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
        self.baseURL = 'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query='
        self.baseDB = 'GPFCore.dbo'
        self.MCATDB = 'GR6Plus7.dbo'
        self.formatURL = '&format=json&timeout={}'

        self.formatURL = '&format=json&timeout={}'

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
        query = (str(self.baseURL)+
            'select vpo.objid, ra, dec, nuv_mag, fuv_mag, fov_radius, nuv_skybg,'
            ' fuv_skybg, nuv_fwhm_world, fuv_fwhm_world, vpe.fexptime,'
            ' vpe.nexptime, fuv_mag_aper_1, fuv_mag_aper_2, fuv_mag_aper_3,'
            ' fuv_mag_aper_4, fuv_mag_aper_5, fuv_mag_aper_6, fuv_mag_aper_7,'
            ' nuv_mag_aper_1, nuv_mag_aper_2, nuv_mag_aper_3, nuv_mag_aper_4,'
            ' nuv_mag_aper_5, nuv_mag_aper_6, nuv_mag_aper_7,'
            ' fuv_magerr_aper_1, fuv_magerr_aper_2, fuv_magerr_aper_3,'
            ' fuv_magerr_aper_4, fuv_magerr_aper_5, fuv_magerr_aper_6,'
            ' fuv_magerr_aper_7, nuv_magerr_aper_1, nuv_magerr_aper_2,'
            ' nuv_magerr_aper_3, nuv_magerr_aper_4, nuv_magerr_aper_5,'
            ' nuv_magerr_aper_6, nuv_magerr_aper_7'
            ' from '+str(self.MCATDB)+'.visitphotoobjall as vpo'
            ' inner join '+str(self.MCATDB)+'.visitphotoextract'
            ' as vpe on vpo.photoextractid=vpe.photoextractid inner join'
            ' '+str(self.MCATDB)+'.fGetNearbyVisitObjEq('+repr(float(self.ra0))+','+
            repr(float(self.dec0))+', '+str(self.radius*60.)+
            ') as nb on vpo.objid=nb.objid'+str(self.formatURL))
        self.assertEqual(gq.exposure_ranges(self.NUV,self.ra0,self.dec0,t0=self.t0,t1=self.t1,detsize=self.detsize),query)

    def test_exposure_range(self):
        self.assertEqual(gq.exposure_range(self.NUV,self.ra0,self.dec0,t0=self.t0,t1=self.t1),"http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select startTimeRange, endTimeRange from fGetTimeRanges(766525332,866526576,176.919525856,0.255696872807) where band='NUV'&format=json&timeout={}")

    def test_aperture(self):
        self.assertEqual(gq.aperture(self.NUV,self.ra0,self.dec0,self.t0,self.t1,self.radius),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=Select  sum(photonCount) from dbo.fGetNearbyObjEqCountNUV(176.919525856,0.255696872807,0.004,766525332995,866526576995,0)&format=json&timeout={}')

    def test_deadtime1(self):
        self.assertEqual(gq.deadtime1(self.NUV,self.t0,self.t1),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select count(*) from NUVPhotonsV where time between 766525332995 and 866526576995&format=json&timeout={}')

    def test_deadtime2(self):
        self.assertEqual(gq.deadtime2(self.NUV,self.t0,self.t1),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select count(*) from NUVPhotonsNULLV where time between 766525332995 and 866526576995&format=json&timeout={}')

    def test_deadtime(self):
        self.assertEqual(gq.deadtime(self.NUV,self.t0,self.t1),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select sum(dt)*0.0000057142857142857145 / (866526576.995-766525332.995) from(select count(*) as dt from NUVPhotonsNULLV where time between 766525332995 and 866526576995 union all select count(*) as dt from NUVPhotonsV where time between 766525332995 and 866526576995) x&format=json&timeout={}')

    def test_boxcount(self):
        self.assertEqual(gq.boxcount(self.NUV,self.t0,self.t1,self.xr,self.yr),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select count(*) from NUVPhotonsNULLV where time between 766525332995 and 866526576995 and x between 200 and 400 and y between 300 and 500&format=json&timeout={}')

    def test_stimcount(self):
        margin=[90.01,90.01]
        avgstim = CalUtils.avg_stimpos(self.NUV,self.eclipse)
        query = ('{baseURL}select count(*) from {baseDB}.{band}PhotonsNULLV '+
            'where time between {t0} and {t1} and ('+
            '((x between {x10} and {x11}) and (y between {y10} and {y11})) or '+
            '((x between {x20} and {x21}) and (y between {y20} and {y21})) or '+
            '((x between {x30} and {x31}) and (y between {y30} and {y31})) or '+
            '((x between {x40} and {x41}) and (y between {y40} and {y41}))'+
            '){formatURL}').format(baseURL=self.baseURL,
                baseDB=self.baseDB, band=self.NUV,
                t0=str(long(t0*self.tscale)), t1=str(long(t1*self.tscale)),
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
        self.assertEqual(gq.boxcentroid(self.NUV,self.t0,self.t1,self.xr,self.yr),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select avg(x), avg(y) from NUVPhotonsNULLV where time between 766525332995 and 866526576995 and x between 200 and 400 and y between 300 and 500&format=json&timeout={}')

    def test_boxtimes(self):
        self.assertEqual(gq.boxtimes(self.NUV,self.t0,self.t1,self.xr,self.yr),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select time from NUVPhotonsNULLV where time between 766525332995 and 866526576995 and x between 200 and 400 and y between 300 and 500&format=json&timeout={}')

    def test_centroid(self):
        self.assertEqual(gq.centroid(self.NUV,self.ra0,self.dec0,self.t0,self.t1,self.radius),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select avg(ra), avg(dec) from NUVPhotonsV where time between 766525332995 and 866526576995 and ra between 176.91552585600002 and 176.923525856 and dec between 0.251696872807 and 0.259696872807&format=json&timeout={}')

    def test_allphotons(self):
        self.assertEqual(gq.allphotons(self.NUV,self.ra0,self.dec0,self.t0,self.t1,self.radius),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select time,ra,dec,xi,eta from dbo.fGetNearbyObjEqNUVAllColumns(176.919525856,0.255696872807,0.004,766525332995,866526576995,0)&format=json&timeout={}')

    def test_shutter(self):
        self.assertEqual(gq.shutter(self.NUV,self.t0,self.t1),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select shutter*0.05 from fGetNUVShutter(766525332995,866526576995)&format=json&timeout={}')

    def test_shutdead(self):
        self.assertEqual(gq.shutdead(self.NUV,self.t0,self.t1),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=SELECT shutter*0.05 FROM fGetNUVShutter(766525332995,866526576995) AS time UNION ALL SELECT SUM(dt) * 0.0000057142857142857145 / (866526576.995-766525332.995) AS dead FROM(SELECT count(*) AS dt FROM NUVPhotonsNULLV WHERE time BETWEEN 766525332995 AND 866526576995 UNION ALL SELECT count(*) AS dt FROM NUVPhotonsV WHERE time BETWEEN 766525332995 AND 866526576995) x&format=json&timeout={}')

    def test_exptime(self):
        self.assertEqual(gq.exptime(self.NUV,self.t0,self.t1),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select * from fGetNUVEffectiveExposureTime(766525332995,866526576995,1.0)&format=json&timeout={}')

    def test_aspect(self):
        self.assertEqual(gq.aspect(self.t0,self.t1),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select eclipse, filename, time, ra, dec, twist, flag, ra0, dec0, twist0 from aspect where time between 766525332995 and 866526576995 order by time&format=json&timeout={}')

    def test_aspect_ecl(self):
        self.assertEqual(gq.aspect_ecl(self.eclipse),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select eclipse, filename, time, ra, dec, twist, flag, ra0, dec0, twist0 from aspect where eclipse=23456 order by time&format=json&timeout={}')

    def test_box(self):
        self.assertEqual(gq.box(self.NUV,self.ra0,self.dec0,self.t0,self.t1,self.radius),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select time,ra,dec from NUVPhotonsV where time between 766525332995 and 866526576995 and ra between 176.91552585600002 and 176.923525856 and dec between 0.251696872807 and 0.259696872807 and flag=0&format=json&timeout={}')

    def test_rect(self):
        self.assertEqual(gq.rect(self.NUV,self.ra0,self.dec0,self.t0,self.t1,self.radius,self.radius),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select time,ra,dec from fGetObjFromRectNUV(176.917525856,176.92152585600002,0.253696872807,0.257696872807,766525332995,866526576995,0)&format=json&timeout={}')

suite = unittest.TestLoader().loadTestsFromTestCase(TestGQueryFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)
