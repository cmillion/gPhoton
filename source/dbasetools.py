# Contains tools for deriving data from the database that are used by
#  a number of different modules.
import numpy as np
import gQuery
from MCUtils import print_inline,area

# Find the contiguous time ranges within a time range at a specific location.
def fGetTimeRanges(band,skypos,trange=[1,1000000000000],tscale=1000.,detsize=1.25,verbose=0,maxgap=1.,minexp=1.):
	try:
		times = np.array(gQuery.getArray(gQuery.exposure_ranges(band,skypos[0],skypos[1],t0=trange[0],t1=trange[1],detsize=detsize,tscale=tscale),verbose=verbose),dtype='float64')[:,0]/tscale
	except:
		return np.array([],dtype='float64')
	if verbose:
		print_inline('Parsing '+str(len(times)-1)+' seconds of exposure.')
	blah = []
	for i in xrange(len(times[0:-1])):
		blah.append(times[i+1]-times[i])
	# A drop in data with duration greater than maxgap initiates a
	#  new exposure
	gaps = np.where(np.array(blah)>maxgap)
	ngaps = len(gaps[0])
	chunks = []
	for i in range(ngaps):
		if not i:
			chunk = [times[0],times[gaps[0][i]]]
		elif i==ngaps-1:
			chunk = [times[gaps[0][i]+1],times[-1]]
		else:
			chunk = [times[gaps[0][i]+1],times[gaps[0][i+1]]]
		# If the duration of this slice is less than minexp, do not
		#  count it as valid exposure.
		if chunk[1]-chunk[0]<minexp:
			continue
		else:
			chunks.append(chunk)
	if not ngaps:
		if times.min()==times.max():
			chunks.append([times.min(),times.min()+1])
		else:
			chunks.append([times.min(),times.max()])

	return np.array(chunks,dtype='float64')

def compute_exptime(band,trange,verbose=0,skypos=[False,False],detsize=1.25):
	if skypos[0] and skypos[1]:
		tranges = fGetTimeRanges(band,skypos,verbose=verbose,trange=trange)
	else:
		tranges = [trange]
	exptime = 0.
	for trange in tranges:
		if trange[0]==trange[1]:
			continue
		rawexpt  = trange[1]-trange[0]
		shutdead = gQuery.getArray(gQuery.shutdead(band,trange[0],trange[1]),verbose=verbose)
		exptime += (rawexpt-shutdead[0][0])*(1.-shutdead[1][0])

	return exptime

# Get the MCAT estimated sky background at a location
def mcat_skybg(band,skypos,radius,verbose=0):
	# Setting maglimit to 30 so that it gets _everything_.
	sources = gQuery.getArray(gQuery.mcat_sources(band,skypos[0],skypos[1],radius,maglimit=30))

	# The MCAT reports skybg in photons/sec/sq.arcsec
	if band=='NUV':
		skybg = np.float32(np.array(sources)[:,5]).mean()
	else:
		skybg = np.float32(np.array(sources)[:,6]).mean()

	# And radius is in degrees
	return skybg*area(radius*60.*60.)

