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
        self.baseURL = 'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query='
        self.formatURL = '&format=json&timeout={}'

    def test_baseURL(self):
        self.assertEqual(gq.baseURL,self.baseURL)

    def test_formatURL(self):
        self.assertEqual(gq.formatURL,self.formatURL)

    def test_mcat_sources(self):
        self.assertEqual(gq.mcat_sources(self.NUV,self.ra0,self.dec0,self.radius,maglimit=self.maglimit),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select ra, dec, nuv_mag, fuv_mag, fov_radius, nuv_skybg, fuv_skybg, nuv_fwhm_world, fuv_fwhm_world from Gr6plus7.Dbo.photoobjall as p inner join Gr6plus7.Dbo.photoextract as pe on p.photoextractid=pe.photoextractid inner join gr6plus7.dbo.fgetnearbyobjeq(176.919525856, 0.255696872807, 1.44) as nb on p.objid=nb.objid and (band=3 or band=1) and NUV_mag<30.0&format=json&timeout={}')

    def test_exposure_ranges(self):
        self.assertEqual(gq.exposure_ranges(self.NUV,self.ra0,self.dec0,t0=self.t0,t1=self.t1,detsize=self.detsize),"http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select distinct time from fGetNearbyAspectEq(176.919525856,0.255696872807,((1.25/2.0)*60.0),766525332995,866526576995) where band='NUV' or band='FUV/NUV' order by time&format=json&timeout={}")

    def test_exposure_range(self):
        self.assertEqual(gq.exposure_range(self.NUV,self.ra0,self.dec0,t0=self.t0,t1=self.t1),"http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select startTimeRange, endTimeRange from fGetTimeRanges(766525332,866526576,176.919525856,0.255696872807) where band='NUV'&format=json&timeout={}")

    def test_aperture(self):
        self.assertEqual(gq.aperture(self.NUV,self.ra0,self.dec0,self.t0,self.t1,self.radius),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=Select  sum(photonCount) from dbo.fGetNearbyObjEqCountNUV(176.919525856,0.255696872807,0.004,766525332995,866526576995,0)&format=json&timeout={}')

    def test_deadtime1(self):
        self.assertEqual(gq.deadtime1(self.NUV,self.t0,self.t1),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select count(*) from NUVPhotonsV where time between 766525332995 and 866526576995&format=json&timeout={}')

    def test_deadtime2(self):
        self.assertEqual(gq.deadtime2(self.NUV,self.t0,self.t1),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select count(*) from NUVPhotonsNULLV where time between 766525332995 and 866526576995&format=json&timeout={}')

    def test_hyper_deadtime(self):
        self.assertEqual(gq.hyper_deadtime(self.NUV,[[self.t0,self.t1]]),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select sum(dt) * 0.0000057142857142857145 / (866526576.995-766525332.995) from(select count(*) as dt from NUVPhotonsNULLV where time between 766525332995 and 866526576995 union all select count(*) as dt from NUVPhotonsV where time between 766525332995 and 866526576995) x&format=json&timeout={}')

    def test_deadtime(self):
        self.assertEqual(gq.deadtime(self.NUV,self.t0,self.t1),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select sum(dt) * 0.0000057142857142857145 / (866526576.995-766525332.995) from(select count(*) as dt from NUVPhotonsNULLV where time between 766525332995 and 866526576995 union all select count(*) as dt from NUVPhotonsV where time between 766525332995 and 866526576995) x&format=json&timeout={}')

    def test_boxcount(self):
        self.assertEqual(gq.boxcount(self.NUV,self.t0,self.t1,self.xr,self.yr),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select count(*) from NUVPhotonsNULLV where time between 766525332995 and 866526576995 and x between 200 and 400 and y between 300 and 500&format=json&timeout={}')

    def test_stimcount(self):
        self.assertEqual(gq.stimcount(self.NUV,self.t0,self.t1,eclipse=self.eclipse),'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select count(*) from NUVPhotonsNULLV where time between 766525332995 and 866526576995 and ((x between -40902.9566054 and -38284.6717091 and y between 34381.2426431 and 36999.5275393) or (x between 34598.6815899 and 37216.9664861 and y between 34379.6427578 and 36997.9276541) or (x between -40897.1388409 and -38278.8539446 and y between -38627.3380359 and -36009.0531396) or (x between 34613.3714451 and 37231.6563414 and y between -38656.7177464 and -36038.4328502))&format=json&timeout={}')

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

