import unittest
import gQuery as gq

class TestgQueryFunctions(unittest.TestCase):
	def setUp(self):
		self.band   = 'NUV'
		self.ra0    = 176.919525856
		self.dec0   = 0.255696872807
		self.t0     = 766525332.995
		self.t1	    = 866526576.995
		self.radius = 0.004
		self.eclipse= 23456
		self.baseURL='http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query='
		self.formatURL='&format=json&timeout={}'

	def test_exposure_range(self):
		self.assertEqual(gq.exposure_range(self.band,self.ra0,self.dec0),self.baseURL+'select startTimeRange, endTimeRange from fGetTimeRanges(1,10000000000000,176.919525856,0.255696872807) where band=\'NUV\''+self.formatURL)

	def test_aperture(self):
		self.assertEqual(gq.aperture(self.band,self.ra0,self.dec0,self.t0,self.t1,self.radius),self.baseURL+'Select  sum(photonCount) from dbo.fGetNearbyObjEqCountNUV(176.919525856,0.255696872807,0.004,766525332995,866526576995,0)'+self.formatURL)

	def test_deadtime1(self):
		self.assertEqual(gq.deadtime1(self.band,self.t0,self.t1),self.baseURL+'select count(*) from NUVPhotons where time between 766525332995 and 866526576995'+self.formatURL)

	def test_deadtime2(self):
		self.assertEqual(gq.deadtime2(self.band,self.t0,self.t1),self.baseURL+'select count(*) from NUVPhotonsNULL where time between 766525332995 and 866526576995'+self.formatURL)

	def test_shutter(self):
		#query = gQuery.shutter(self.band,self.t0,self.t1)
		self.assertTrue(True)

	def test_aspect(self):
		self.assertEqual(gq.aspect(self.t0,self.t1),self.baseURL+'select eclipse, filename, time, ra, dec, twist, flag, ra0, dec0, twist0 from aspect where between 766525332995 and 866526576995 order by time'+self.formatURL)

	def test_aspect_ecl(self):
		self.assertEqual(gq.aspect_ecl(self.eclipse),self.baseURL+'select eclipse, filename, time, ra, dec, twist, flag, ra0, dec0, twist0 from aspect where eclipse=23456 order by time'+self.formatURL)

	def test_box(self):
		self.assertEqual(gq.box(self.band,self.ra0,self.dec0,self.t0,self.t1,self.radius),self.baseURL+'select time,ra,dec from NUVPhotons where time between 766525332995 and 866526576995 and ra between 176.915525856 and 176.923525856 and dec between 0.251696872807 and 0.259696872807 and flag=0'+self.formatURL)

	def test_make_int(self):
		self.assertTrue(1==1)

#       def tearDown():
#

suite = unittest.TestLoader().loadTestsFromTestCase(TestgQueryFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)

