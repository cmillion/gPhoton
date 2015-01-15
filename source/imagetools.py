import gQuery
import numpy as np
import MCUtils as mc
from astropy import wcs as pywcs
from astropy.io import fits as pyfits
import scipy.misc
import scipy.special # erfc
import scipy.ndimage
from FileUtils import flat_filename
import gnomonic
import dbasetools as dbt
import galextools as gxt

def define_wcs(skypos,skyrange,width=False,height=False,verbose=0,
			   pixsz=0.000416666666666667):
	"""Define the world coordinate system (WCS)."""
	if verbose:
		mc.print_inline('Defining World Coordinate System (WCS).')
	wcs = pywcs.WCS(naxis=2) # NAXIS = 2
	imsz = gxt.deg2pix(skypos,skyrange)
	wcs.wcs.cdelt = np.array([-pixsz,pixsz])
	wcs.wcs.ctype = ['RA---TAN','DEC--TAN']
	wcs.wcs.crpix = [(imsz[1]/2.)+0.5,(imsz[0]/2.)+0.5]
	wcs.wcs.crval = skypos
	return wcs

def movie_tbl(band,tranges,verbose=0,framesz=0,retries=20):
	"""Initialize a FITS table to contain movie frame information."""
	if verbose:
		mc.print_inline('Populating exposure time table.')
	tstarts,tstops,exptimes=[],[],[]
	for trange in tranges:
		stepsz = framesz if framesz else trange[1]-trange[0]
		steps = np.ceil((trange[1]-trange[0])/stepsz)
		for i,t0 in enumerate(np.arange(trange[0],trange[1],stepsz)):
			t1 = trange[1] if i==steps else t0+stepsz
			tstarts.append(t0)
			tstops.append(t1)
			exptimes.append(dbt.compute_exptime(band,[t0,t1],
							verbose=verbose,retries=retries))
	col1 = pyfits.Column(name='tstart',format='E',array=np.array(tstarts))
	col2 = pyfits.Column(name='tstop',format='E',array=np.array(tstops))
	col3 = pyfits.Column(name='exptime',format='E',array=np.array(exptimes))
	cols = pyfits.ColDefs([col1,col2,col3])
	tbl  = pyfits.BinTableHDU.from_columns(cols)

	return tbl

def fits_header(band,skypos,tranges,skyrange,width=False,height=False,
				verbose=0,tscale=1000.,hdu=False,retries=20):
	"""Populate a FITS header."""
	if verbose:
		mc.print_inline('Populating FITS header.')
	hdu = hdu if hdu else pyfits.PrimaryHDU()
	wcs = define_wcs(skypos,skyrange,width=width,height=height)
	hdu.header['CDELT1'],hdu.header['CDELT2'] = wcs.wcs.cdelt
	hdu.header['CTYPE1'],hdu.header['CTYPE2'] = wcs.wcs.ctype
	hdu.header['CRPIX1'],hdu.header['CRPIX2'] = wcs.wcs.crpix
	hdu.header['CRVAL1'],hdu.header['CRVAL2'] = wcs.wcs.crval
	#hdu.header['RA_CENT'],hdu.header['DEC_CENT'] = wcs.wcs.crval # Dupe.
	hdu.header['EQUINOX'],hdu.header['EPOCH'] = 2000., 2000.
	hdu.header['BAND'] = 1 if band=='NUV' else 2
	# Do we want to set the following?
	#hdu.header['OW'] = 1
	#hdu.header['DIRECT'] = 1
	#hdu.header['GRISM'] = 0
	#hdu.header['OPAQUE'] = 0

	# Put the total exposure time into the primary header
	hdu.header['EXPTIME'] = 0.
	for trange in tranges:
		hdu.header['EXPTIME'] += dbt.compute_exptime(band,trange,
												verbose=verbose,retries=retries)

	if len(tranges)==1:
	# Put the time range into the primary header for a single frame image
		hdu.header['EXPSTART'],hdu.header['EXPEND'] = tranges[0]
		# These are the proper keywords for this:
		hdu.header['TIME-OBS'],hdu.header['TIME-END'] = tranges[0]

	return hdu

def sigmaclip_bg(data,radius,annulus,skypos,maxiter=10,sigmaclip=3.,
				 gausslim=50.,verbose=0,pixsz=0.000416666666666667):
	"""Produce an estimate of background counts within an aperture (radius)
	using a sigma clipping method for extracting the background from an
	annulus.
	This attempts to reproduce the calcuations of the backcalc() function in
	mapaps/poissonbg.c of the mission pipeline. (Probably written by Ted Wyder.)
	"""
	# FIXME: Does not apply response!
	d = mc.angularSeparation(
		skypos[0],skypos[1],data['ra'],data['dec'])
	# This cut is now handled by the ugly loop below, which barely dodges a
	# conceptula issue about fractional pixels...
	#ix = np.where((d>annulus[0]) & (d<annulus[1]))
	coo=zip(data['ra'],data['dec'])

	imsz=gxt.deg2pix(skypos,[annulus[1]*2,annulus[1]*2])
	wcs=define_wcs(skypos,[annulus[1]*2,annulus[1]*2])
	foc=wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
	H,xedges,yedges=np.histogram2d(foc[:,1]-0.5,foc[:,0]-0.5,bins=imsz,
										range=([ [0,imsz[0]],[0,imsz[1]] ]))

	# Convert Gaussian sigma to a probability
	problim = 0.5*scipy.special.erfc(sigmaclip/np.sqrt(2.0));

	# Mask out non-annulus regions... there's probalby a more pythonic way
	bgimg=np.copy(H)
	for i in range(H.shape[0]):
		for j in range(H.shape[1]):
			# Add a little buffer to account for pixel widths?
			# FIXME? including everything within the annulus...
#			if (mc.distance(H.shape[0]/2.,H.shape[1]/2.,i,j)<annulus[0]/pixsz or
			if	mc.distance(H.shape[0]/2.,H.shape[1]/2.,i,j)>annulus[1]/pixsz:#):

				bgimg[i,j]=-1

	ix=np.where(bgimg>=0)
	m,s=bgimg[ix].mean(),bgimg[ix].std()
	d = 1.
	for i in range(maxiter):
		if d<=10e-5 or m<2:
			continue
		if m>=gausslim:
			# Mask anything outside of 3 sigma from the mean (of unmasked data)
			klim=m+sigmaclip*np.sqrt(m)#s
			klo=m-sigmaclip*np.sqrt(m)#s
			if verbose:
				print 'Gaussian cut: {klo} to {klim}'.format(klo=klo,klim=klim)
		else:
			klim = scipy.special.gammainccinv(m,problim)
			klo = -1 # None
			if verbose:
				print 'Poisson cut: {klo} to {klim}'.format(klo=klo,klim=klim)
		ix = np.where((bgimg>=klim) | (bgimg<=klo))
		bgimg[ix]=-1
		ix=np.where(bgimg>=0)
		d = np.abs((bgimg[ix].mean()-m)/m)# - 1)
		m,s=bgimg[ix].mean(),bgimg[ix].std()
	ix = np.where(bgimg>=0)
	return mc.area(radius)*bgimg[ix].mean()/mc.area(pixsz)

def countmap(band,skypos,tranges,skyrange,width=False,height=False,
			 verbose=0,tscale=1000.,memlight=False,hdu=False,retries=20):
	"""Create a count (cnt) map."""
	imsz = gxt.deg2pix(skypos,skyrange)
	count = np.zeros(imsz)
	for trange in tranges:
		# If memlight is requested, break the integration into
		#  smaller chunks.
		step = memlight if memlight else trange[1]-trange[0]
		for i in np.arange(trange[0],trange[1],step):
			t0,t1=i,i+step
			if verbose:
				mc.print_inline('Coadding '+str(t0)+' to '+str(t1))
			events = gQuery.getArray(gQuery.rect(band,skypos[0],skypos[1],t0,t1,
												 skyrange[0],skyrange[1]),
									 verbose=verbose,retries=retries)

			# Check that there is actually data here.
			if not events:
				if verbose>1:
					print "No data in "+str([t0,t1])
				continue

			times = np.array(events,dtype='float64')[:,0 ]/tscale
			coo   =	np.array(events,dtype='float64')[:,1:]

			# If there's no data, return a blank image.
			if len(coo)==0:
				if verbose:
					print 'No data in this frame: '+str([t0,t1])
				continue

			# Define World Coordinate System (WCS)
			wcs = define_wcs(skypos,skyrange,width=False,height=False)

			# Map the sky coordinates onto the focal plane
			foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)

			# Bin the events into actual image pixels
			H,xedges,yedges=np.histogram2d(foc[:,1]-0.5,foc[:,0]-0.5,
								bins=imsz,range=([ [0,imsz[0]],[0,imsz[1]] ]))
			count += H

	return count

def write_jpeg(filename,band,skypos,tranges,skyrange,width=False,height=False,
			   stepsz=1.,clobber=False,verbose=0,tscale=1000.,retries=20):
	"""Write a 'preview' jpeg image from a count map."""
	scipy.misc.imsave(filename,countmap(band,skypos,tranges,skyrange,
					  width=width,height=height,verbose=verbose,tscale=tscale,
					  retries=retries))
	return

def rrhr(band,skypos,tranges,skyrange,width=False,height=False,stepsz=1.,
		 verbose=0,calpath='../cal/',tscale=1000.,response=True,hdu=False,
		 retries=20):
	"""Generate a high resolution relative response (rrhr) map."""
	imsz = gxt.deg2pix(skypos,skyrange)
	# TODO the if width / height

	flat = mc.get_fits_data(flat_filename(band,calpath),verbose=verbose)
	flatinfo = mc.get_fits_header(flat_filename(band,calpath))
	npixx,npixy 	= flat.shape
	fltsz 		= flat.shape
	pixsz = flatinfo['CDELT2']
	detsize = 1.25

	# Rotate the flat into the correct orientation to start.
	flat   = np.flipud(np.rot90(flat))

	# NOTE: This upsample interpolation is done _last_ in the canonical
	#	pipeline as part of the poissonbg.c routine.
	# 	The interpolation function is "congrid" in the same file.
	# TODO: Should this be first order interpolation? (i.e. bilinear)
	hrflat = scipy.ndimage.interpolation.zoom(flat,4.,order=0,prefilter=False)
	img = np.zeros(hrflat.shape)[
				hrflat.shape[0]/2.-imsz[0]/2.:hrflat.shape[0]/2.+imsz[0]/2.,
				hrflat.shape[1]/2.-imsz[1]/2.:hrflat.shape[1]/2+imsz[1]/2.]

	for trange in tranges:
		t0,t1=trange
		entries = gQuery.getArray(gQuery.aspect(t0,t1),retries=retries)
		n = len(entries)

		asptime = np.float64(np.array(entries)[:,2])/tscale
		aspra   = np.float32(np.array(entries)[:,3])
		aspdec  = np.float32(np.array(entries)[:,4])
		asptwist= np.float32(np.array(entries)[:,5])
		aspflags= np.float32(np.array(entries)[:,6])
		asptwist= np.float32(np.array(entries)[:,9])
		aspra0  = np.zeros(n)+skypos[0]
		aspdec0 = np.zeros(n)+skypos[1]

		xi_vec, eta_vec = gnomonic.gnomfwd_simple(
							aspra,aspdec,aspra0,aspdec0,-asptwist,1.0/36000.,0.)

		col = 4.*( ((( xi_vec/36000.)/(detsize/2.)*(detsize/(fltsz[0]*pixsz)) + 1.)/2. * fltsz[0]) - (fltsz[0]/2.) )
		row = 4.*( (((eta_vec/36000.)/(detsize/2.)*(detsize/(fltsz[1]*pixsz)) + 1.)/2. * fltsz[1]) - (fltsz[1]/2.) )

		vectors = mc.rotvec(np.array([col,row]),-asptwist)

		for i in range(n):
			if verbose>1:
				mc.print_inline('Stamping '+str(asptime[i]))
				# FIXME: Clean this mess up a little just for clarity.
	        	img += scipy.ndimage.interpolation.shift(scipy.ndimage.interpolation.rotate(hrflat,-asptwist[i],reshape=False,order=0,prefilter=False),[vectors[1,i],vectors[0,i]],order=0,prefilter=False)[hrflat.shape[0]/2.-imsz[0]/2.:hrflat.shape[0]/2.+imsz[0]/2.,hrflat.shape[1]/2.-imsz[1]/2.:hrflat.shape[1]/2+imsz[1]/2.]*dbt.compute_exptime(band,[asptime[i],asptime[i]+1],verbose=verbose,retries=retries)*gxt.compute_flat_scale(asptime[i]+0.5,band,verbose=0)

	return img

def movie(band,skypos,tranges,skyrange,framesz=0,width=False,height=False,
		  verbose=0,tscale=1000.,memlight=False,coadd=False,
		  response=False,calpath='../cal/',hdu=False,retries=20):
	"""Generate a movie (mov)."""
	# Not defining stepsz effectively creates a count map.
	mv = []
	rr = []
	if coadd:
		if verbose>2:
			print 'Coadding across '+str(tranges)
		mv.append(countmap(band,skypos,tranges,skyrange,width=width,
				  height=height,verbose=verbose,tscale=tscale,memlight=memlight,
				  hdu=hdu,retries=retries))
		rr.append(rrhr(band,skypos,tranges,skyrange,response=response,width=width,height=height,stepsz=1.,verbose=verbose,calpath=calpath,tscale=tscale,hdu=hdu,retries=retries)) if response else rr.append(np.ones(np.shape(mv)[1:]))
	else:
		for trange in tranges:
			stepsz = framesz if framesz else trange[1]-trange[0]
			steps = np.ceil((trange[1]-trange[0])/stepsz)
			for i,t0 in enumerate(np.arange(trange[0],trange[1],stepsz)):
				if verbose>1:
					mc.print_inline('Movie frame '+str(i+1)+' of '+
																str(int(steps)))
				t1 = trange[1] if i==steps else t0+stepsz
				mv.append(countmap(band,skypos,[[t0,t1]],skyrange,width=width,
								   height=height,verbose=verbose,tscale=tscale,
								   memlight=memlight,hdu=hdu,retries=retries))
	# FIXME: This should not create an rr unless it's requested...
				rr.append(rrhr(band,skypos,[[t0,t1]],skyrange,response=response,width=width,height=height,stepsz=1.,verbose=verbose,calpath=calpath,tscale=tscale,retries=retries)) if response else rr.append(np.ones(np.shape(mv)[1:]))

	return np.array(mv),np.array(rr)

def create_images(band,skypos,tranges,skyrange,framesz=0,width=False,
				  height=False,verbose=0,tscale=1000.,memlight=False,
				  coadd=False,response=False,calpath='../cal/',hdu=False,
				  retries=20):
	count,rr = movie(band,skypos,tranges,skyrange,framesz=framesz,
					 width=width,height=height,verbose=verbose,tscale=tscale,
					 memlight=memlight,coadd=coadd,response=response,
					 calpath=calpath,hdu=hdu,retries=retries)
	intensity = []
	for i in range(count.shape[0]):
		int_temp = count[i]/rr[i] # Might divide by zero in rare cases...
		cut = np.where(np.isfinite(int_temp))
		int_clean = np.zeros(int_temp.shape)
		int_clean[cut]=int_temp[cut]
		intensity.append(int_clean.tolist())

	return np.array(count),np.array(rr),np.array(intensity)

def write_images(band,skypos,tranges,skyrange,write_cnt=False,write_int=False,
				 write_rr=False,framesz=0,width=False,height=False,verbose=0,
				 tscale=1000.,memlight=False,coadd=False,response=False,
				 calpath='../cal/',clobber=False,retries=20):
	"""Generate a write various maps to files."""
	# No files were requested, so don't bother doing anything.
	if not (write_cnt or write_int or write_rr):
		return
	count,rr,intensity = create_images(band,skypos,tranges,skyrange,
									   framesz=framesz,width=width,
									   height=height,verbose=verbose,
									   tscale=tscale,memlight=memlight,
									   coadd=coadd,response=response,
									   calpath=calpath,retries=retries)
	# Add a conditional so that this is only created for multi-frame images
	tbl = movie_tbl(band,tranges,framesz=framesz,verbose=verbose,
															retries=retries)
	if write_cnt:
		hdu = pyfits.PrimaryHDU(count)
		hdu = fits_header(band,skypos,tranges,skyrange,width=width,
						  height=height,verbose=verbose,tscale=tscale,hdu=hdu,
						  retries=retries)
		hdulist = pyfits.HDUList([hdu,tbl])
		if verbose:
			print 'Writing count image to '+str(write_cnt)
		hdulist.writeto(write_cnt,clobber=clobber)
	if write_rr:
		hdu = pyfits.PrimaryHDU(rr)
		hdu = fits_header(band,skypos,tranges,skyrange,width=width,
						  height=height,verbose=verbose,tscale=tscale,
						  hdu=hdu,retries=retries)
		hdulist = pyfits.HDUList([hdu,tbl])
		if verbose:
			print 'Writing response image to '+str(write_rr)
                hdulist.writeto(write_rr,clobber=clobber)
	if write_int:
		hdu = pyfits.PrimaryHDU(intensity)
		hdu = fits_header(band,skypos,tranges,skyrange,width=width,
						  height=height,verbose=verbose,tscale=tscale,hdu=hdu,
						  retries=retries)
		hdulist = pyfits.HDUList([hdu,tbl])
		if verbose:
			print 'Writing intensity image to '+str(write_int)
		hdulist.writeto(write_int,clobber=clobber)

	return
