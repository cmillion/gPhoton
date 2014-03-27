import CalUtils
import unittest
#import unittest2
#from CalibrationTools import load_txy
import numpy as np
import os

class TestCalUtilsFunctions(unittest.TestCase):

	def setUp(self):
		self.NUV = 'NUV'
		self.FUV = 'FUV'
		# definse an eclipse that is before the CSP
		self.preCSP  = 31000
		# defines an eclipse that is after the CSP
		self.postCSP = 40000
		self.detsize = 1.5
		self.imsz    = 3840.
		self.pixsz   = 0.00016666667
		self.npixx   = 3840.
		self.npixy   = 3840.
		self.fill    = self.detsize/(self.npixx*self.pixsz)
		self.raw6 = '../e31000/MISWZN11_12494_0315_0002-fd-raw6.fits'
		self.ssd  = 'testSSD.tbl'
		self.calpath = '../cal/'
		self.pltscl = 68.754932
		self.aspum  = self.pltscl/1000.

	def test_clk_cen_scl_slp(self):
		self.assertEqual(CalUtils.clk_cen_scl_slp(self.NUV,self.preCSP), (2007.0, 1992.0, 7400.0, 6070.0, 8.79, 14.01, 0.53, 0.0))
                self.assertEqual(CalUtils.clk_cen_scl_slp(self.FUV,self.preCSP), (1997.0, 1993.0, 7200.0, 6670.0, 7.78, 10.73, 0.0, 0.0))
                self.assertEqual(CalUtils.clk_cen_scl_slp(self.NUV,self.postCSP),(2007.0, 2016.0, 7400.0, 6070.0, 8.79, 14.01, 0.53, 0.0))
                self.assertEqual(CalUtils.clk_cen_scl_slp(self.FUV,self.postCSP),(1997.0, 1993.0, 7200.0, 6670.0, 7.78, 10.73, 0.0, 0.0))


	def test_xieta2colrow(self):
		col,row = CalUtils.xieta2colrow(0.,0.,self.detsize,self.fill,self.npixx,self.npixy)
		self.assertAlmostEqual(col,1920.)
		self.assertAlmostEqual(row,1920.)

	def test_avg_stimpos(self):
		self.assertEqual(CalUtils.avg_stimpos(self.FUV,self.preCSP), {'x1':-2541.88,'x2':2632.06,'x3':-2541.53,'x4':2631.68,'y1':2455.28,'y2':2455.02,'y3':-2550.89,'y4':-2550.92})
                self.assertEqual(CalUtils.avg_stimpos(self.FUV,self.postCSP),{'x1':-2541.88,'x2':2632.06,'x3':-2541.53,'x4':2631.68,'y1':2455.28,'y2':2455.02,'y3':-2550.89,'y4':-2550.92})
		self.assertEqual(CalUtils.avg_stimpos(self.NUV,self.preCSP), {'x1':-2722.27,'x2':2468.84,'x3':-2721.87,'x4':2469.85,'y1':2453.89,'y2':2453.78,'y3':-2565.81,'y4':-2567.83})
		self.assertEqual(CalUtils.avg_stimpos(self.NUV,self.postCSP),{'x1':-2722.53,'x2':2470.29,'x3':-2721.98,'x4':2471.09,'y1':2549.96,'y2':2550.10,'y3':-2538.57,'y4':-2538.62})

	def test_find_stims_index(self):
		t_sham  = np.sort(np.random.rand(self.npixx*self.npixx*4/100)*10+918992428.362)
		# Spread photon events evenly over detector space
		detpts = range(-int(self.npixx*10),int(self.npixx*10),100)
		x_sham = []
		y_sham = []
		for i in detpts:
			for j in detpts:
				x_sham.append(i)
				y_sham.append(j)

		# FUV
		index1,index2,index3,index4 = CalUtils.find_stims_index(x_sham,y_sham,self.FUV,self.preCSP,margin=90.001)
		self.assertTrue(len(index1)==676)
		self.assertTrue(len(index2)==378)
		self.assertTrue(len(index3)==702)
		self.assertTrue(len(index4)==378)

		# NUV, preCSP
                index1,index2,index3,index4 = CalUtils.find_stims_index(x_sham,y_sham,self.NUV,self.preCSP,margin=90.001)
                self.assertTrue(len(index1)==52)
                self.assertTrue(len(index2)==702)
                self.assertTrue(len(index3)==48)
                self.assertTrue(len(index4)==624)

		# NUV, postCSP
                index1,index2,index3,index4 = CalUtils.find_stims_index(x_sham,y_sham,self.NUV,self.postCSP,margin=90.001)
                self.assertTrue(len(index1)==52)
                self.assertTrue(len(index2)==676)
                self.assertTrue(len(index3)==52)
                self.assertTrue(len(index4)==676)

	def test_find_stims(self):
		t_sham  = np.sort(np.random.rand(self.npixx*self.npixx*4/100)*10+918992428.362)
		# Spread photon events evenly over detector space
		detpts = range(-int(self.npixx*10),int(self.npixx*10),100)
		x_sham = []
		y_sham = []
		for i in detpts:
			for j in detpts:
				x_sham.append(i)
				y_sham.append(j)

		# FUV
		t,x,y,ix = CalUtils.find_stims(t_sham,x_sham,y_sham,self.FUV,self.preCSP)
		self.assertTrue(len(t)==2134)
		self.assertTrue(len(x)==len(t) and len(y)==len(t) and len(ix)==len(t))

		# NUV, preCSP
		t,x,y,ix = CalUtils.find_stims(t_sham,x_sham,y_sham,self.NUV,self.preCSP)
                self.assertTrue(len(t)==1426)
                self.assertTrue(len(x)==len(t) and len(y)==len(t) and len(ix)==len(t))

		# NUV, postCSP
                t,x,y,ix = CalUtils.find_stims(t_sham,x_sham,y_sham,self.NUV,self.postCSP)
                self.assertTrue(len(t)==1456)
                self.assertTrue(len(x)==len(t) and len(y)==len(t) and len(ix)==len(t))

        def test_create_SSD(self):
                CalUtils.create_SSD(self.raw6,self.FUV,self.preCSP,ssdfile=self.ssd)
		self.assertTrue(os.path.exists(self.ssd))
                os.remove(self.ssd)
		self.assertFalse(os.path.exists(self.ssd))

	def test_get_stim_coefs(self):
		# Generate a test SSD file
                CalUtils.create_SSD(self.raw6,self.FUV,self.preCSP,ssdfile=self.ssd)

		coef0,coef1 = CalUtils.get_stim_coefs(self.ssd)
		self.assertTrue(abs(coef0+449716.93841278262)<15000.)
		self.assertTrue(abs(coef1-0.00049489726822987549)*10000<1.)

		# Remove the test SSD file
                os.remove(self.ssd)

	# DEPRECATED
	#def test_compute_aspect_vectors(self):
	#	p = [0.,0.,0.,0.]
	#	v = np.array([0.,0.,0.])
	#	e = CalUtils.compute_aspect_vectors(p,p,p,p,p,p,p)
	#	self.assertTrue(e[0][0]==0. and e[0][1]==0. and e[0][2]==0.)
        #       self.assertTrue(e[1][0]==0. and e[1][1]==0. and e[1][2]==0.)
	#	p = [30.,30.,30.,30.]
	#	e = CalUtils.compute_aspect_vectors(p,p,p,p,p,p,p)
        #       self.assertTrue(e[0][0]==0. and e[0][1]==0. and e[0][2]==0.)
        #       self.assertTrue(e[1][0]==0. and e[1][1]==0. and e[1][2]==0.)

	def test_find_FUV_offset(self):
		self.assertTrue(1==1)

	def test_postCSP_caldata(self):
		wig2, wig2data, wlk2, wlk2data, clk2, clk2data = CalUtils.postCSP_caldata(self.calpath)
		self.assertTrue(wig2.size==262144)
		self.assertTrue(wig2data=={'start':-2000,'inc':40})
		self.assertTrue(wlk2.size==25600)
		self.assertTrue(wlk2data=={'start':-2000,'inc':40})
		self.assertTrue(clk2.size==800)
		self.assertTrue(clk2data=={'start':-2000,'inc':40})

	def test_rtaph_yap(self):
		ya = np.array([0])
		yb,yamc=ya,ya
		self.assertEqual(CalUtils.rtaph_yap(ya,yb,yamc)[0],0)
		ya = np.array([4])
		yb,yamc=ya,ya
		self.assertEqual(CalUtils.rtaph_yap(ya,yb,yamc)[0],4)

	def test_rtaph_yac(self):
		self.assertTrue(1==1)

	def test_rtaph_yac2(self):
		wig2, wig2data, wlk2, wlk2data, clk2, clk2data = CalUtils.postCSP_caldata(self.calpath)
		# some random test data pulled from e31000
		q  = np.array([ 15,       10,    10,     16     ])
		xb = np.array([ 5,        4,     4,      3      ])
		yb = np.array([ 2,        4,     4,      3      ])
		ya = np.array([ 6,        9,     9,      27     ])
		y  = np.array([ -15622.88,8905.9,7350.05,1030.08])

		yac2 = CalUtils.rtaph_yac2(q,xb,yb,ya,y,self.calpath,self.aspum,wig2, wig2data, wlk2, wlk2data, clk2, clk2data)
		self.assertAlmostEqual(yac2[0], 19.74403814)
		self.assertAlmostEqual(yac2[1],-30.22328638)
		self.assertAlmostEqual(yac2[2],-32.91836577)
		self.assertAlmostEqual(yac2[3], 3.28412804 )

	def test_raw6_to_stims(self):
		stim1,stim2,stim3,stim4=CalUtils.raw6_to_stims(self.raw6,self.FUV,self.preCSP)
		self.assertTrue(len(stim1['t'])==27253)
		self.assertTrue(len(stim2['t'])==27063)
		self.assertTrue(len(stim3['t'])==27308)
		self.assertTrue(len(stim4['t'])==27232)

	def test_compute_stimstats(self):
		Mx,Bx,My,By,stimsep,yac=CalUtils.compute_stimstats(self.raw6,self.FUV,self.preCSP)
		self.assertAlmostEqual(Mx,  0.99997150415456593)
		self.assertAlmostEqual(Bx,  0.55909143454953347)
		self.assertAlmostEqual(My,  0.99999376556508601)
		self.assertAlmostEqual(By, -1.0915058186972255 )
		self.assertAlmostEqual(stimsep, 5089.9043198566624)
		self.assertTrue(len(yac)==40)

#	def tearDown():
#

suite = unittest.TestLoader().loadTestsFromTestCase(TestCalUtilsFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)

