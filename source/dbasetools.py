# Contains tools for working with data from the database that are used by
#  a number of different modules.
import numpy as np
import gQuery
from MCUtils import print_inline,area,distance,angularSeparation
from galextools import GPSSECS, zpmag

def get_aspect(band,skypos,trange=[6e8,11e8],tscale=1000.,verbose=0,
               detsize=1.25):
    """Get aspect solution in a dict() for given time range."""
    asp = np.array(gQuery.getArray(gQuery.aspect_skypos(skypos[0],skypos[1],
                   detsize=detsize),verbose=verbose))
    data = {'eclipse':np.array(asp[:,0],dtype='int16'),'filename':asp[:,1],
            't':np.array(asp[:,2],dtype='float64')/tscale,
            'ra':np.array(asp[:,3],dtype='float64'),
            'dec':np.array(asp[:,4],dtype='float64'),
            'twist':np.array(asp[:,5],dtype='float64'),
            'flag':np.array(asp[:,6],dtype='int8'),
            'ra0':np.array(asp[:,7],dtype='float64'),
            'dec0':np.array(asp[:,8],dtype='float64'),
            'twist0':np.array(asp[:,9],dtype='float64')}
    ix = np.where((data['t']>trange[0]) & (data['t']<trange[1]) &
                  (angularSeparation(skypos[0],skypos[1],
                                     data['ra'],data['dec'])<=detsize/2.))
    for key in data.keys():
        data[key] = data[key][ix]
    return data

def fGetTimeRanges(band,skypos,trange=None,tscale=1000.,detsize=1.25,verbose=0,
                   maxgap=1.,minexp=1.,retries=100.,predicted=False):
    """Find the contiguous time ranges within a time range at a
    specific location.
	minexp - Do not include exposure time less than this.
	maxgap - Gaps in exposure longer than this initiate a new time range.
	detsize - Fiddle with this if you want to exlude the edges of the
    detector.
    predicted - Use the aspect solutions to estimate what exposure will be
    available once the database is fully populated.
	"""
    try:
        if not np.array(trange).tolist():
            trange = [1,1000000000000]
        if len(np.shape(trange))==2:
            trange=trange[0]
        if verbose and predicted:
            print "Querying coverage of GR6/7. (Not current database.)"
        times = (np.array(gQuery.getArray(gQuery.exposure_ranges(band,
            skypos[0],skypos[1],t0=trange[0],t1=trange[1],detsize=detsize,
            tscale=tscale),verbose=verbose,retries=retries),
                                                dtype='float64')[:,0]/tscale
            if not predicted else get_aspect(band,skypos,trange,detsize=detsize,
                                        tscale=tscale,verbose=verbose)['t'])
    except IndexError:
        if verbose:
            print "No exposure time available at {pos}".format(pos=skypos)
        return np.array([],dtype='float64')
    except TypeError:
        print "Is one of the inputs malformed?"
        raise
    except:
        raise
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

def exposure(band,trange,verbose=0,retries=20,photonfile=None):
    """Compute the effective exposure time for a time range."""
    rawexpt = trange[1]-trange[0]
    if rawexpt<=0:
        return 0.
    if photonfile:
        print ""
    else:
        shutter = gQuery.getArray(gQuery.shutdead(band,trange[0],trange[1]),
                                        verbose=verbose,retries=retries)[0][0]
    # NOTE: The deadtime correction in shutdead does not work properly in FUV
    # so we're doing it separately for now.
    deadtime = gQuery.getValue(gQuery.deadtime(band,trange[0],trange[1]),
                                            verbose=verbose,retries=retries)
    #return (rawexpt-shutdead[0][0])*(1.-shutdead[1][0])
    return (rawexpt-shutter)*(1.-deadtime)

def compute_exptime(band,trange,verbose=0,skypos=None,detsize=1.25,
                    retries=20,chunksz=10.e6,coadd=False):
    """Compute the effective exposure time."""
    # FIXME: This skypos[] check appears to not work properly and leads
    #  to dramatic _underestimates_ of the exposure time.
    if skypos and coadd:
        tranges = fGetTimeRanges(band,skypos,verbose=verbose,trange=trange,
                                 retries=retries)
    else:
        tranges=[trange]
    if verbose>1:
        print_inline('Computing exposure within {tr}'.format(tr=tranges))
        if skypos and coadd:
            print 'Based on skypos {sp}'.format(sp=skypos)
    exptime = 0.
    for trange in tranges:
        # Traverse the exposure ranges in manageable chunks to hopefully keep
        # the query response time below the HTTP timeout...
        chunks = (np.linspace(trange[0],trange[1],
                             num=np.ceil((trange[1]-trange[0])/chunksz)) if
                                 (trange[1]-trange[0])>chunksz else
                                 np.array(trange))
        for i,t in enumerate(chunks[:-1]):
            exptime += exposure(band,[chunks[i],chunks[i+1]],verbose=verbose,
                                retries=retries)
    return exptime

def get_mcat_data(skypos,rad):
    out = np.array(gQuery.getArray(
                         gQuery.mcat_visit_sources(skypos[0],skypos[1],rad)))
    # FIXME: The APER entries should really be generated
    try:
        return {'objid':np.array(out[:,0],dtype='int64'),
            'ra':np.array(out[:,1],dtype='float32'),
            'dec':np.array(out[:,2],dtype='float32'),
            'NUV':{'mag':np.array(out[:,3],dtype='float32'),
                   'skybg':np.array(out[:,6],dtype='float32'),
                   'expt':np.array(out[:,11],dtype='float32'),
                   'fwhm':np.array(out[:,8],dtype='float32'),
                   1:{'mag':np.array(out[:,19],dtype='float32')+zpmag('NUV'),
                      'err':np.array(out[:,33],dtype='float32')},
                   2:{'mag':np.array(out[:,20],dtype='float32')+zpmag('NUV'),
                      'err':np.array(out[:,34],dtype='float32')},
                   3:{'mag':np.array(out[:,21],dtype='float32')+zpmag('NUV'),
                      'err':np.array(out[:,35],dtype='float32')},
                   4:{'mag':np.array(out[:,22],dtype='float32')+zpmag('NUV'),
                      'err':np.array(out[:,36],dtype='float32')},
                   5:{'mag':np.array(out[:,23],dtype='float32')+zpmag('NUV'),
                      'err':np.array(out[:,37],dtype='float32')},
                   6:{'mag':np.array(out[:,24],dtype='float32')+zpmag('NUV'),
                      'err':np.array(out[:,38],dtype='float32')},
                   7:{'mag':np.array(out[:,25],dtype='float32')+zpmag('NUV'),
                      'err':np.array(out[:,39],dtype='float32')} },
            'FUV':{'mag':np.array(out[:,4],dtype='float32'),
                   'skybg':np.array(out[:,7],dtype='float32'),
                   'expt':np.array(out[:,10],dtype='float32'),
                   'fwhm':np.array(out[:,9],dtype='float32'),
                   1:{'mag':np.array(out[:,12],dtype='float32')+zpmag('FUV'),
                      'err':np.array(out[:,26],dtype='float32')},
                   2:{'mag':np.array(out[:,13],dtype='float32')+zpmag('FUV'),
                      'err':np.array(out[:,27],dtype='float32')},
                   3:{'mag':np.array(out[:,14],dtype='float32')+zpmag('FUV'),
                      'err':np.array(out[:,28],dtype='float32')},
                   4:{'mag':np.array(out[:,15],dtype='float32')+zpmag('FUV'),
                      'err':np.array(out[:,29],dtype='float32')},
                   5:{'mag':np.array(out[:,16],dtype='float32')+zpmag('FUV'),
                      'err':np.array(out[:,30],dtype='float32')},
                   6:{'mag':np.array(out[:,17],dtype='float32')+zpmag('FUV'),
                      'err':np.array(out[:,31],dtype='float32')},
                   7:{'mag':np.array(out[:,18],dtype='float32')+zpmag('FUV'),
                      'err':np.array(out[:,32],dtype='float32')} } }
    except IndexError:
        return False
    except:
        raise

def mcat_skybg(band,skypos,radius,verbose=0,retries=20):
	"""Estimate the sky background using the MCAT skybg for nearby sources."""
	# Setting maglimit to 30 so that it gets _everything_...
	sources = gQuery.getArray(gQuery.mcat_sources(band,skypos[0],skypos[1],
                                                  radius,maglimit=30),
                                                  retries=retries)

	# The MCAT reports skybg in photons/sec/sq.arcsec
	if band=='NUV':
		skybg = np.float32(np.array(sources)[:,5]).mean()
	else:
		skybg = np.float32(np.array(sources)[:,6]).mean()

	# And radius is in degrees
	return skybg*area(radius*60.*60.)

def get_mags(band,ra0,dec0,radius,maglimit,mode='coadd',
                   zpmag={'NUV':20.08,'FUV':18.82},verbose=0):
    """Given RA, Dec and search radius, searches the coadd MCAT for sources.
    Returns a dict() which contains magnitudes for all of the APER settings.
    Note: Visit mode returns a lot more sources, more slowly than coadd mode
    given the same search parameters. You should probably use smaller search
    radii in visit mode. If you're just trying to find unique sources in a
    large region, use coadd mode and then pass the result through the
    parse_unique_sources() function contained in this module.
    """
    zpf,zpn = zpmag['FUV'],zpmag['NUV']
    if mode=='coadd':
        out =np.array(gQuery.getArray(
                  gQuery.mcat_sources(band,ra0,dec0,radius,maglimit=maglimit),
                  verbose=verbose))
        if not len(out):
            print "Warning: No sources found!"
            return 0
        return {'ra':out[:,0],'dec':out[:,1],
                'FUV':{'mag':out[:,3],1:out[:,9]+zpf,2:out[:,10]+zpf,
                       3:out[:,11]+zpf,4:out[:,12]+zpf,5:out[:,13]+zpf,
                       6:out[:,14],7:out[:,15]+zpf},
                'NUV':{'mag':out[:,2],1:out[:,16]+zpn,2:out[:,17]+zpn,
                       3:out[:,18]+zpn,4:out[:,19]+zpn,5:out[:,20]+zpn,
                       6:out[:,21]+zpn,7:out[:,22]+zpn}}
    elif mode=='visit':
        out = np.array(gQuery.getArray(
                       gQuery.mcat_visit_sources(ra0,dec0,radius),
                       verbose=verbose))
        # NOTE: For runtime considerations, mcat_visit_sources() does not
        # make any slices on time or maglimit, so we need to do it here.
        ix = np.where((out[:,2 if band=='NUV' else 3]<maglimit) &
                                           (out[:,2 if band=='NUV' else 3]>0))
        return {'ra':out[:,0][ix],'dec':out[:,1][ix],
                'NUV':{'mag':out[:,2][ix],'expt':out[:,8][ix],
                1:out[:,18][ix]+zpn,2:out[:,19][ix]+zpn,3:out[:,20][ix]+zpn,
                4:out[:,21][ix]+zpn,5:out[:,22][ix]+zpn,6:out[:,23][ix]+zpn,
                7:out[:,24][ix]+zpn},
                'FUV':{'mag':out[:,3][ix],'expt':out[:,9][ix],
                1:out[:,11][ix]+zpf,2:out[:,12][ix]+zpf,3:out[:,13][ix]+zpf,
                4:out[:,14][ix]+zpf,5:out[:,15][ix]+zpf,6:out[:,16][ix]+zpf,
                7:out[:,17][ix]+zpf}}
    else:
        print "mode must be in [coadd,visit]"
        return None

def exp_from_objid(objid):
    out = np.array(gQuery.getArray(gQuery.mcat_objid_search(objid)))
    return {'NUV':{'expt':np.array(out[:,7],dtype='float')[0],'t0':np.array(out[:,9],dtype='float64')[0]-GPSSECS,'t1':np.array(out[:,10],dtype='float64')[0]-GPSSECS},'FUV':{'expt':np.array(out[:,8],dtype='float')[0],'t0':np.array(out[:,11],dtype='float64')[0]-GPSSECS,'t1':np.array(out[:,12],dtype='float64')[0]-GPSSECS}}

def parse_unique_sources(ras,decs,fmags,nmags,margin=0.005):
    """Iteratively returns unique sources based upon a _margin_ within
    which two sources should be considered the same sources. Is a little
    bit sensitive to the first entry and could probably be written to be
    more robust, but works well enough.
    """
    skypos = zip(ras,decs)
    for i,pos in enumerate(skypos):
        ix = np.where(angularSeparation(pos[0],pos[1],ras,decs)<=margin)
        skypos[i]=[ras[ix].mean(),decs[ix].mean()]
        a = skypos #unique_sources(data['ra'],data['dec'])
    b = []
    for i in a:
        if not (i in b):
            b+=[i]
    return b

def find_unique_sources(band,ra0,dec0,searchradius,maglimit=20.0,margin=0.001,
                                                                    verbose=0):
    coadds = get_mags(band,ra0,dec0,searchradius,maglimit,mode='coadd',
                                                            verbose=verbose)
    if not coadds:
        return None
    else:
        return np.array(parse_unique_sources(coadds['ra'],coadds['dec'],
                    coadds['FUV']['mag'],coadds['NUV']['mag'],margin=margin))

def avg_sources(band,skypos,radius=0.001,maglimit=20.0,verbose=0,
                                                    catalog='MCAT',retries=20):
	"""Return the mean position of sources within the search radius."""
	out = np.array(gQuery.getArray(gQuery.mcat_sources(band,skypos[0],
        skypos[1],radius,maglimit=maglimit),verbose=verbose,retries=retries))
	ix = np.where(out[:,-2]>0) if band=='NUV' else np.where(out[:,-1]>0)
	fwhm = out[ix,-2].mean() if band=='NUV' else out[ix,-1].mean()
	return out[ix,0].mean(),out[ix,1].mean(),round(fwhm,4)

def nearest_source(band,skypos,radius=0.01,maglimit=20.0,verbose=0,
                                                    catalog='MCAT',retries=20):
    """Return targeting parameters for the nearest MCAT source to a position.
    """
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

    #dist=np.sqrt( (out[:,0]-skypos[0])**2 + (out[:,1]-skypos[1])**2)
    dist = angularSeparation(out[:,0],out[:,1],skypos[0],skypos[1])
    if verbose > 1:
        print "Finding nearest among "+str(len(dist))+" nearby sources."
    # Note that this doesn't cope with multiple entries for the same source.
    s = out[np.where(dist == dist.min())][0]
    # RA, Dec, NUV mag, FUV mag, NUV fwhm, FUV fwhm
    return avg_sources(band,[s[0],s[1]],verbose=verbose,retries=retries)
    #return s[0],s[1],s[2],s[3],s[7],s[8]

def nearest_distinct_source(band,skypos,radius=0.1,maglimit=20.0,verbose=0,
                                                    catalog='MCAT',retries=20):
    """Return parameters for the nearest non-targeted source."""
    out = np.array(gQuery.getArray(gQuery.mcat_sources(band,skypos[0],skypos[1],radius,maglimit=maglimit),verbose=verbose,retries=retries))
    #dist=np.sqrt( (out[:,0]-skypos[0])**2 + (out[:,1]-skypos[1])**2)
    dist = angularSeparation(out[:,0],out[:,1],skypos[0],skypos[1])
    ix = np.where(dist>0.005)
    return np.array(out)[ix][np.where(dist[ix]==dist[ix].min())][0]

def suggest_bg_radius(band,skypos,radius=0.1,maglimit=20.0,verbose=0,
                                                    catalog='MCAT',retries=20):
    """Returns a recommended background radius based upon the
    positions and FWHM of nearby sources in the MCAT.
    """
    nearest = nearest_distinct_source(band,skypos,verbose=verbose,retries=retries)
    #dist = np.sqrt( (nearest[0]-skypos[0])**2 + (nearest[1]-skypos[1])**2 )
    dist = angularSeparation(nearest[0],nearest[1],skypos[0],skypos[1])
    return round(dist-3*nearest[-2 if band=='NUV' else -1],4)

def optimize_annulus(optrad,outann,verbose=0):
	"""Suggest optiumum annulus dimensions."""
	if outann<=round(2*optrad,4):
		print "Warning: There are known sources within the background annulus."
		print "Use --hrbg to mask these out. (Will increase run times.)"
	return round(1.2*optrad,4),round(2*optrad,4)

def suggest_parameters(band,skypos,verbose=0,retries=20):
    mcat = get_mcat_data(skypos,0.0005)
    ix = np.where((mcat[band]['mag']>0) & (mcat[band]['fwhm']>0))
    pos,fwhm = None, None
    if mcat['objid'].any(): # There is a known star at the target position!
        pos = [mcat['ra'][ix].mean(),mcat['dec'][ix].mean()]
        radius = 2*mcat[band]['fwhm'][ix].mean()
        if verbose:
            print 'Recentering on {pos}.'.format(pos=pos)
            print 'Using aperture radius of {rad} degrees.'.format(rad=fwhm)
    else: # There is no known star at the target position...
        pos = skypos
        radius = gt.aper2deg(4)
    annulus = [3*radius, 5*radius]

    return pos[0],pos[1],radius,annulus[0],annulus[1]
