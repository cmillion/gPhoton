# Computing deadtime:
#  Two main ways to estimate the deadtime, each of which is implemented in the
#  compute_deadtime() funciton in CalibrationTools.py. The first "correct" way
#  compares the measured stim count rate against the nominal commanded rate.
#  This doesn't work very well for short time scales, though, because it's
#  subject to Poisson statistics. The second "less correct" way uses a conversion
#  fact between global detector count rate and deadtime. The second method is
#  far more consistent across exposure times (because the total number of
#  detector counts, even per second, is quite high.
#
#  The implementation in compute_deadtime() works straight from the photon CSV
#  files so we will have to adapt them to work with the database.
#
#  It would also be useful to recompute the scale fact for Method 2. Patrick
#  Morissey, who did the original work with Tim Conrow, said that the scale
#  factor is slightly different in each band, but they forced them to be
#  the same for consistency. They also only used a subset of the total mission
#  data, and he doesn't have any of their original code left.
#
#  We should test the possibility that the scale factor for the global dead time
#  may have changed over the mission.

from CalUtils import avg_stimpos
from dbasetools import fGetTimeRanges, compute_exptime
import numpy as np
import gQuery
from MCUtils import print_inline

band = 'NUV'
t0,t1 = 766525332.995,766526576.995
skypos = [176.919525856024,0.255696872807351]

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

skypos = [294.86493, 47.30736]
tranges = fGetTimeRanges(band,skypos)
deadtime = [netdead(band,trange[0],trange[1]) for trange in tranges]

skypos = [323.06766667,0.25400000]
tranges = fGetTimeRanges(band,skypos)
deadtime0 = [netdead(band,trange[0],trange[1]) for trange in tranges[np.random.random_integers(0,tranges.shape[0]-1,100)]]

%pylab
# Two methods vs. global count rate
plt.plot(np.array(deadtime)[:,-1]/(np.array(deadtime)[:,1]-np.array(deadtime)[:,0]),np.array(deadtime)[:,2],'.')
plt.plot(np.array(deadtime)[:,-1]/(np.array(deadtime)[:,1]-np.array(deadtime)[:,0]),np.array(deadtime)[:,3],'.')

# Two methods vs. each other
plt.plot(np.array(deadtime)[:,2],np.array(deadtime)[:,3],'.')

