# Contains tools for deriving data from the database that are used by
#  a number of different modules.
import numpy as np
import gQuery
from MCUtils import print_inline,area

def fGetTimeRanges(band,skypos,trange=None,tscale=1000.,detsize=1.25,verbose=0,maxgap=1.,minexp=1.,retries=100.):
    """Find the contiguous time ranges within a time range at a 
    specific location.
	minexp - Do not include exposure time less than this.
	maxgap - Gaps in exposure longer than this initiate a new time range.
	detsize - Fiddle with this if you want to exlude the edges of the
    detector.
	"""
    try:
        #FIXME: t[01] appears to have no impact on this
        # if trange is not set, set it to an arbitrary large range in order
        # to capture the whole mission
        if not trange:
            trange = [1,1000000000000]
        if len(np.shape(trange))==2:
            trange=trange[0]
        times = np.array(gQuery.getArray(gQuery.exposure_ranges(band,skypos[0],skypos[1],t0=trange[0],t1=trange[1],detsize=detsize,tscale=tscale),verbose=verbose,retries=retries),dtype='float64')[:,0]/tscale
    except:
        return np.array([],dtype='float64')
    if verbose:
        print_inline('Parsing '+str(len(times)-1)+' seconds of exposure.: ['+str(trange[0])+', '+str(trange[1])+']')
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

def exposure(band,trange,verbose=0,retries=20):
    """Compute the effective exposure time for a time range."""
    rawexpt = trange[1]-trange[0]
    if rawexpt<=0:
        return 0.
    shutdead = gQuery.getArray(gQuery.shutdead(band,trange[0],trange[1]),verbose=verbose,retries=retries)
    return (rawexpt-shutdead[0][0])*(1.-shutdead[1][0])

def compute_exptime(band,trange,verbose=0,skypos=None,detsize=1.25,retries=20):
    """Compute the effective exposure time."""
    # FIXME: This skypos[] check appears to not work properly and leads
    #  to dramatic _underestimates_ of the exposure time.
    if skypos:
        tranges = fGetTimeRanges(band,skypos,verbose=verbose,trange=trange,
                                 retries=retries)
    else:
        tranges=[trange]
    exptime = 0.
    for trange in tranges:
        # To create manageable queries, only compute exposure time in
        # chunks <=10e6 seconds
        chunksz = 10.e6
        chunks = (np.linspace(trange[0],trange[1],
                             num=np.ceil((trange[1]-trange[0])/chunksz)) if 
                                 (trange[1]-trange[0])>chunksz else
                                 np.array(trange))
        for i,t in enumerate(chunks[:-1]):
            exptime += exposure(band,[chunks[i],chunks[i+1]],verbose=verbose,
                                retries=retries)
    return exptime

def mcat_skybg(band,skypos,radius,verbose=0,retries=20):
	"""Estimate the sky background using the MCAT skybg for nearby sources."""
	# Setting maglimit to 30 so that it gets _everything_.
	sources = gQuery.getArray(gQuery.mcat_sources(band,skypos[0],skypos[1],radius,maglimit=30),retries=retries)

	# The MCAT reports skybg in photons/sec/sq.arcsec
	if band=='NUV':
		skybg = np.float32(np.array(sources)[:,5]).mean()
	else:
		skybg = np.float32(np.array(sources)[:,6]).mean()

	# And radius is in degrees
	return skybg*area(radius*60.*60.)

def avg_sources(band,skypos,radius=0.001,maglimit=22.0,verbose=0,catalog='MCAT',retries=20):
	"""Return the mean position of sources within the search radius."""
	out = np.array(gQuery.getArray(gQuery.mcat_sources(band,skypos[0],skypos[1],radius,maglimit=maglimit),verbose=verbose,retries=retries))
	ix = np.where(out[:,-2]>0) if band=='NUV' else np.where(out[:,-1]>0)
	fwhm = out[ix,-2].mean() if band=='NUV' else out[ix,-1].mean()
	return out[ix,0].mean(),out[ix,1].mean(),round(fwhm,4)

def nearest_source(band,skypos,radius=0.01,maglimit=22.0,verbose=0,catalog='MCAT',retries=20):
	"""Return targeting parameters for the nearest MCAT source to a position."""
	out = np.array(gQuery.getArray(gQuery.mcat_sources(band,skypos[0],skypos[1],radius,maglimit=maglimit),verbose=verbose,retries=retries))
	if not len(out) and band=='FUV':
		if verbose:
			print "No nearby MCAT source found in FUV. Trying NUV..."
		band = 'NUV'
		out = np.array(gQuery.getArray(gQuery.mcat_sources(band,skypos[0],skypos[1],radius,maglimit=maglimit),verbose=verbose,retries=retries))
	if not len(out) and band=='NUV':
		if verbose:
			print "No nearby MCAT source found. Using input sky position."
		return skypos[0],skypos[1],0.01
	
	dist=np.sqrt( (out[:,0]-skypos[0])**2 + (out[:,1]-skypos[1])**2)
	if verbose > 1:
			print "Finding nearest among "+str(len(dist))+" nearby sources."
	# Note that this doesn't cope with multiple entries for the same source.
	s = out[np.where(dist == dist.min())][0]
	# RA, Dec, NUV mag, FUV mag, NUV fwhm, FUV fwhm
	return avg_sources(band,[s[0],s[1]],verbose=verbose,retries=retries)
	#return s[0],s[1],s[2],s[3],s[7],s[8]

def nearest_distinct_source(band,skypos,radius=0.1,maglimit=22.0,verbose=0,catalog='MCAT',retries=20):
	"""Return parameters for the nearest non-targeted source."""
	out = np.array(gQuery.getArray(gQuery.mcat_sources(band,skypos[0],skypos[1],radius,maglimit=maglimit),verbose=verbose,retries=retries))
	dist=np.sqrt( (out[:,0]-skypos[0])**2 + (out[:,1]-skypos[1])**2)
	ix = np.where(dist>0.005)
	return np.array(out)[ix][np.where(dist[ix]==dist[ix].min())][0]

def suggest_bg_radius(band,skypos,radius=0.1,maglimit=22,verbose=0,catalog='MCAT',retries=20):
	"""Returns a recommended background radius based upon the
	positions and FWHM of nearby sources in the MCAT.
	"""
	nearest = nearest_distinct_source(band,skypos,verbose=verbose,retries=retries)
	dist = np.sqrt( (nearest[0]-skypos[0])**2 + (nearest[1]-skypos[1])**2 )
	return round(dist-3*nearest[-2 if band=='NUV' else -1],4)

def optimize_annulus(optrad,outann,verbose=0):
	"""Suggest optiumum annulus dimensions."""
	if outann<=round(2*optrad,4):
		print "Warning: There are known sources within the background annulus."
		print "Use --hrbg to mask these out. (Will increase run times.)"
#		if verbose:
#			print "Warning: There is no optimum background annulus."
#			print "Using no background correction; CAVEAT EMPTOR!"
#		inann,outann=0,0
#	elif outann<=1.1*optrad:
#		if verbose:
#			print "Warning: Background of questionable utility."
#		inann=optrad
#	elif outann>=3*optrad:
#		inann,outann=round(1.5*optrad,4),round(3*optrad,4)
#	elif outann < 3*optrad:
#		inann=round(1.1*optrad,4)
#	else:
#		print "I DON'T THINK THIS SHOULD HAPPEN!!&!*!!!"
#	return inann,outann
	# Doing the fancy thing above resulted in too many [0,0] annuli
	return round(1.2*optrad,4),round(2*optrad,4)

def suggest_parameters(band,skypos,radius=0.01,maglimit=22.0,verbose=0,catalog='MCAT',retries=20):
	"""Suggest an optimum aperture position and size."""
	ra,dec,fwhm=nearest_source(band,skypos,radius=radius,maglimit=maglimit,verbose=verbose,retries=retries)
	optrad = round(2*fwhm,4)
	outann = suggest_bg_radius(band,skypos,maglimit=maglimit,verbose=verbose,retries=retries)
	annulus = optimize_annulus(optrad,outann,verbose=verbose)
	if verbose:
		print "Suggested sky position [RA,Dec]: ["+str(ra)+", "+str(dec)+"]"
		print "Suggested aperture radius (deg): "+str(optrad)
		print "Suggested background annulus:    ["+str(annulus[0])+", "+str(annulus[1])+"]"
	return ra,dec,optrad,annulus[0],annulus[1]

