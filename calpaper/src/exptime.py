from CalUtils import avg_stimpos
from dbasetools import fGetTimeRanges, compute_exptime
import numpy as np
import gQuery
from MCUtils import print_inline

def netdead(band,t0,t1,tstep=1.,refrate=79.,verbose=0):
	print 'Time range: ['+str(t0)+', '+str(t1)+']'
	refrate     = 79. # counts per second, nominal stim rate
	feeclkratio = 0.966 # not sure what detector property this adjusts for
	tec2fdead   = 5.52e-6 # TEC to deadtime correction (Method 2)
	stimcount = gQuery.getValue(gQuery.stimcount(band,t0,t1))
	totcount = (gQuery.getValue(gQuery.deadtime1(band,t0,t1)) +
		    gQuery.getValue(gQuery.deadtime2(band,t0,t1)))
	exptime = t1-t0
	# Method 0
	# Empirical formula
	dead0 = tec2fdead*(totcount/exptime)/feeclkratio
	minrate,maxrate = refrate*.4,refrate+2.
	# Method 2
	# Direct measurement of stims
	bins = np.linspace(0.,exptime-exptime%tstep,exptime//tstep+1)+t0
	h = np.zeros(len(bins))
	for i,t in enumerate(bins):
		print_inline(t1-t)
		h[i] = gQuery.getValue(gQuery.stimcount(band,t,t+tstep-0.0001))
	
	#h,xh = np.histogram(stimt-trange[0],bins=bins)
	ix = ((h<=maxrate) & (h>=minrate)).nonzero()[0]
	dead2 = (1.-((h[ix]/tstep)/feeclkratio)/refrate).mean()
	# Method 1
	# Empirical formula except with counts binned into 1s
	h = np.zeros(len(bins))
	for i,t in enumerate(bins):
		print_inline(t1-t)
		h[i] = gQuery.getValue(gQuery.deadtime1(band,t,t+tstep-0.0001)) + gQuery.getValue(gQuery.deadtime2(band,t,t+tstep-0.0001))

	dead1 = (tec2fdead*(h/tstep)/feeclkratio).mean()
	expt = compute_exptime(band,[t0,t1])
	print dead0,dead1,dead2
	return [t0, t1, dead0, dead1, dead2, expt, stimcount, totcount]


