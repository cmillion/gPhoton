import csv
import gQuery
from MCUtils import *
from gnomonic import *
from FileUtils import flat_filename
from imagetools import backgroundmap
import dbasetools as dbt
import galextools as gxt

def compute_background(band,skypos,trange,radius,annulus,verbose=0):
    """Estimated counts within an aperture based upon rate within an annulus.
    counts in aperture = (area of aperture) * (counts in annulus) / (area of annulus)
    """
    if annulus[1]-annulus[0]<=0:
        return 0.
    else:
        return (area(radius)/(area(annulus[1])-area(annulus[0]))) * (gQuery.getValue(gQuery.aperture(band,skypos[0],skypos[1],trange[0],trange[1],annulus[1]),verbose=verbose)-gQuery.getValue(gQuery.aperture(band,skypos[0],skypos[1],trange[0],trange[1],annulus[0]),verbose=verbose)) if (annulus[0] and annulus[1]) else 0.

def mask_background_image(band,skypos,trange,radius,annulus,verbose=0,
        maglimit=28.,pixsz=0.000416666666666667,NoData=-999):
    """Creates an image of the sky within an annulus."""
    skyrange = [2*annulus[1], 2*annulus[1]]
    imsz = gxt.deg2pix(skypos,skyrange)
    # FIXME: Duplicate code from backgroundmap()
    # Create a background map with stars blacked out
    bg=backgroundmap(band,skypos,trange,skyrange,verbose=verbose,
                     maglimit=maglimit)
    # Create a position reference array
    xind = np.array([range(int(imsz[1]))]*int(imsz[0])) - (imsz[0]/2.) + 0.5
    yind = np.rot90(
           np.array([range(int(imsz[0]))]*int(imsz[1])) - (imsz[1]/2.)) + 0.5
    distarray = np.sqrt(((xind)**2.) + ((yind)**2.))
    # Cut out the annulus
    # FIXME: (?) Distort the annulus to match the projection?
    ix = np.where(distarray<=annulus[0]/pixsz)
    bg[ix] = NoData
    ix = np.where(distarray>=annulus[1]/pixsz)
    bg[ix] = NoData
    return bg

def compute_background_improved(band,skypos,trange,radius,annulus,verbose=0,
                                maglimit=28.,pixsz=0.000416666666666667,
                                NoData=-999):
    """Estimated counts within the aperture based on counts within the annulus
    with sources within the annulus masked out (to the specified maglimit).
    """
    bg = mask_background_image(band,skypos,trange,radius,annulus,
                               verbose=verbose,maglimit=maglimit,
                               pixsz=pixsz,NoData=NoData)
    # The mean, unmasked background pixel value normalized by the area
    #  of a pixel and scaled to the area of the aperture
    return (area(radius)*bg[np.where(bg>=0)].mean()/(pixsz)**2.)

# This computes the mean response within the aperture. It is slightly less
#  precise but much faster than building a relative response map.
def aperture_response(band,skypos,tranges,radius,verbose=0,tscale=1000.,
                      calpath='../cal/'):
    """Estimates the mean response within the aperture by averaging over
    the portion of the flat upon which the detector falls.
    This might be less accurate than building a response image,
    but it is fast(er).
    """
    if verbose>1:
        print 'Computing relative response over a '+str(radius)+' degree aperture.'
    flat = np.flipud(np.rot90(get_fits_data(flat_filename(band,calpath),
                              verbose=0)))
    flatinfo = get_fits_header(flat_filename(band,calpath))
    detsize, imsz, pixsz = 1.25, flat.shape[0], flatinfo['CDELT2']

    tranges = map(lambda trange: dbt.fGetTimeRanges(band,skypos,verbose=verbose,trange=trange),tranges)[0]

    response = 0.
    # FIXME: Dupe code. A lot of this should be its own function.
    for trange in tranges:
        # Load all the aspect data.
        entries = gQuery.getArray(gQuery.aspect(trange[0],trange[1]),
                                  verbose=verbose)
        n = len(entries)
        asptime	= np.float64(np.array(entries)[:,2])/tscale
        aspra	= np.float32(np.array(entries)[:,3])
        aspdec	= np.float32(np.array(entries)[:,4])
        asptwist= np.float32(np.array(entries)[:,5])
        aspflags= np.float32(np.array(entries)[:,6])
        asptwist= np.float32(np.array(entries)[:,9])
        aspra0	= np.zeros(n)+skypos[0]
        aspdec0	= np.zeros(n)+skypos[1]

        # Compute vectors
        xi_vec, eta_vec = gnomfwd_simple(aspra0,aspdec0,aspra,aspdec,-asptwist,1.0/36000.,0.)
        col  = (((( xi_vec/36000.)/(detsize/2.)*(detsize/(imsz*pixsz))+1.)/2.*imsz)-(imsz/2.))
        row  = ((((eta_vec/36000.)/(detsize/2.)*(detsize/(imsz*pixsz))+1.)/2.*imsz)-(imsz/2.))

        # Build reference arrays in which each element is equal to
        #  its x or y value.
        xind =          np.array([range(int(imsz))]*int(imsz))-(imsz/2.)
        yind = np.rot90(np.array([range(int(imsz))]*int(imsz))-(imsz/2.))

        vectors = rotvec(np.array([col,row]),-asptwist) 
        scales  = gxt.compute_flat_scale(asptime,band,verbose=0)

        pixrad = radius/pixsz#gxt.deg2pix(radius,CDELT2=pixsz)
        # FIXME: How to avoid this loop??!?
        for i in xrange(n):
            if verbose>1:
                print_inline(' Stamping '+str(i+1)+' of '+str(n)+'...')
            asprange=[asptime[i], asptime[i]+1
                      if asptime[i]+1<trange[1] else trange[1]]
            if asprange[0]==asprange[1]:
                if verbose>1:
                    print_inline(str(asprange[0])+'=='+str(asprange[1])+' so skipping...')
                continue
            distarray = np.sqrt(((vectors[0,i]-xind)**2.)+
                                ((vectors[1,i]-yind)**2.))
            ix = np.where(distarray<=pixrad)
            response += dbt.compute_exptime(band,asprange,verbose=verbose)*scales[i]*flat[ix].mean()
            if verbose >1:
                print_inline('                                             ')

    return response

# Returns the flux data within a radius/annulus
def compute_flux(band,skypos,tranges,radius,annulus=[None,None],userr=False,
		 usehrbg=False,verbose=0,calpath='../cal/',detsize=1.25):
	"""Returns a data structure with a bunch of information about the target.
	This includes flux, exposure time, background, response, etc.
	"""
	# Find the exposure time first. If there's not any, don't bother
	#  computing the flux at all.
	data = {'expt':0.,'bg':0., 'bghr':0.,'counts':0.,'rr':1.}
	# Do we really need this? It would be cleaner if it didn't loop.
	for trange in tranges:
		expt = dbt.compute_exptime(band,trange,verbose=verbose,skypos=skypos,detsize=detsize)
		expt = dbt.compute_exptime(band,trange,verbose=verbose,detsize=detsize)
		if not expt:
			continue
		data['expt'] += expt
		data['bg'] += compute_background_improved(band,skypos,trange,radius,annulus,verbose=verbose) if usehrbg else compute_background(band,skypos,trange,radius,annulus,verbose=verbose)
		data['counts'] += gQuery.getValue(gQuery.aperture(band,skypos[0],skypos[1],trange[0],trange[1],radius),verbose=verbose)

	if verbose>1:
		print 'Integrated '+str(data['expt'])+' seconds of exposure.'

	# FIXME: Need to account for 1.018 response shift after 881881215.995
	data['rr'] = aperture_response(band,skypos,tranges,radius,verbose=verbose,calpath=calpath)/data['expt'] if userr and expt>0 else 1

	# TODO: This list of assertions should maybe be folded into a
	#  list comprehension of some kind.
	data['band'] = band
	data['skypos'] = skypos
	data['tranges'] = tranges
	data['tmin'] = np.array(tranges).min()
	data['tmax'] = np.array(tranges).max()
	data['radius'] = radius
	data['annulus'] = annulus
	data['apcor1'] = gxt.apcorrect1(radius,band)
	data['apcor2'] = gxt.apcorrect2(radius,band)
	data['cps'] = ((data['counts']-data['bg'])/data['expt'])/data['rr'] if expt else None
	data['error'] = (np.sqrt(data['counts']-data['bg'])/data['expt'])/data['rr'] if data['cps']>=0 else None
	data['flux'] = gxt.counts2flux(data['cps'],band) if data['cps']>=0 else None
	data['fluxerror'] = gxt.counts2flux(data['error'],band) if data['cps']>=0 else None
	data['mag'] = gxt.counts2mag(data['cps'],band) if data['cps']>0 else None
	if data['cps']>0:
		data['magerr'] = [gxt.counts2mag(data['cps']-data['error'],band) if data['cps']>data['error'] else None, gxt.counts2mag(data['cps']+data['error'],band)]
	else:
		data['magerr']=[None,None]

	return data

def coadd_mag(band,skypos,radius,annulus=[None,None],userr=False,usehrbg=False,
	      verbose=0,detsize=1.25,maxgap=1,minexp=1,calpath='../cal/'):
	tranges = dbt.fGetTimeRanges(band,skypos,verbose=verbose,maxgap=maxgap,minexp=minexp)
	return compute_flux(band,skypos,tranges,radius,annulus=annulus,userr=userr,usehrbg=usehrbg,verbose=verbose,calpath=calpath,detsize=detsize)

def colnames():
	"""Defines the column names for the light curve CSV file."""
	return ['tmin','tmax','radius','exptime','cps','error','flux','fluxerror','mag','magerr','inner annulus','outer annulus','background','response','counts','apcor1','apcor2']

def flux_row_constructor(data):
	"""Constructs rows for the light curve CSV file."""
	return [data['tmin'],data['tmax'],data['radius'],data['expt'],
		data['cps'],data['error'],data['flux'],data['fluxerror'],
		data['mag'],data['magerr'][1],data['annulus'][0],
		data['annulus'][1],data['bg'],data['rr'],data['counts'],
		data['apcor1'],data['apcor2']]

def compute_curve(band,skypos,tranges,radius,annulus=[None,None],stepsz=None,
		  coadd=False,userr=False,usehrbg=False,verbose=0,
		  calpath='../cal/',detsize=1.25,maxgap=1,minexp=1):
	"""Computes a light curve over the requested parameters."""
	data = []
	# If coadd is specified, integrate across all time ranges.
	if coadd:
		# TODO: reduce()
		result = compute_flux(band,skypos,tranges,radius,
				      annulus=annulus,userr=userr,
				      usehrbg=usehrbg,verbose=verbose,
				      calpath=calpath,detsize=detsize)
		data.append(flux_row_constructor(result))
	else:
		# TODO: map()
		for trange in chunks(tranges,length=stepsz,verbose=verbose):
			if verbose>1:
				print 'Computing flux for '+str(trange[0])+' to '+str(trange[1])
			result=compute_flux(band,skypos,[trange],radius,
					    annulus=annulus,verbose=verbose,
					    userr=userr,usehrbg=usehrbg,
					    calpath=calpath)
			if verbose>1:
				print '        '+str(result['mag'])+' magnitude'
                        data.append(flux_row_constructor(result))

	return data

def write_curve(band,skypos,tranges,radius,outfile=False,annulus=[None,None],
		stepsz=None,coadd=False,userr=False,usehrbg=False,verbose=0,
		iocode='wb',calpath='../cal/',detsize=1.25):
	"""Generates a light curve and writes it to the a CSV file."""
	data = compute_curve(band,skypos,tranges,radius,annulus=annulus,
			     stepsz=stepsz,userr=userr,usehrbg=usehrbg,
			     verbose=verbose,coadd=coadd,calpath=calpath,detsize=detsize)
	if outfile:
		spreadsheet = csv.writer(open(outfile,iocode), delimiter=',',
					 quotechar='|',
					 quoting=csv.QUOTE_MINIMAL)
		spreadsheet.writerows(data)
	else:
		for frame in data:
			print frame
	return

# FIXME: This is probably redundant or should be replaced with pandas.
def read_curve(csvfile):
	"""Reads a light curve into a dict."""
	data = {'tstart':np.array([]),'tstop':np.array([]),
		'radius':np.array([]),'exptime':np.array([]),
		'cps':np.array([]),'error':np.array([]),
		'magnitude':np.array([]),'mag_error':np.array([]),
		'inner annulus':np.array([]),'outer annulus':np.array([]),
		'background':np.array([]),'response':np.array([]),
		'counts':np.array([]),'aperture correction 1':np.array([]),
		'aperture correction 2':np.array([])}

	reader = csv.reader(open(csvfile,'rb'),delimiter=',',quotechar='|')
	for row in reader:
		if row[4] == '':
			continue
		data['tstart'] = np.append(data['tstart'],np.float(row[0]))
		data['tstop'] = np.append(data['tstop'],np.float(row[1]))
		data['radius'] = np.append(data['radius'],np.float(row[2]))
		data['exptime'] = np.append(data['exptime'],np.float(row[3]))
		data['cps'] = np.append(data['cps'],np.float(row[4]))
		data['error'] = np.append(data['error'],np.float(row[5]))
		data['magnitude'] = np.append(data['magnitude'],
				    np.float(row[8]))
		data['mag_error'] = np.append(data['mag_error'],
				    np.float(row[9]))
	        #data['inner annulus'] = np.append(data['inner annulus'],np.float(row[10]))
	        #data['outer annulus'] = np.append(data['outer annulus'],np.float(row[11]))
		data['background'] = np.append(data['background'],
				     np.float(row[12]))
		data['response'] = np.append(data['response'],np.float(row[13]))
		data['counts'] = np.append(data['counts'],np.float(row[14]))
		data['aperture correction 1'] = np.append(data['aperture correction 1'],np.float(row[15]))
		data['aperture correction 2'] = np.append(data['aperture correction 2'],np.float(row[16]))

	return data

