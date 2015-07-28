import gQuery
import numpy as np
import MCUtils as mc
from astropy import wcs as pywcs
from astropy.io import fits as pyfits
import scipy.misc
import scipy.special # erfc
import scipy.ndimage
import gnomonic
import dbasetools as dbt
import galextools as gxt
import curvetools as ct
import cal
from gQuery import tscale
from gPhoton import __version__

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
				verbose=0,hdu=False,retries=20):
	"""Populate a FITS header."""
	if verbose:
		mc.print_inline('Populating FITS header.')
	hdu = hdu if hdu else pyfits.PrimaryHDU()
	wcs = define_wcs(skypos,skyrange,width=width,height=height)
	hdu.header['CDELT1'],hdu.header['CDELT2'] = wcs.wcs.cdelt
	hdu.header['CTYPE1'],hdu.header['CTYPE2'] = wcs.wcs.ctype
	hdu.header['CRPIX1'],hdu.header['CRPIX2'] = wcs.wcs.crpix
	hdu.header['CRVAL1'],hdu.header['CRVAL2'] = wcs.wcs.crval
	hdu.header['EQUINOX'],hdu.header['EPOCH'] = 2000., 2000.
	hdu.header['BAND'] = 1 if band=='NUV' else 2
	hdu.header['VERSION'] = 'v{v}'.format(v=__version__)

	# If requested, put the total exposure time into the primary header
	#hdu.header['EXPTIME'] = 0.
	#for trange in tranges:
	#	hdu.header['EXPTIME'] += dbt.compute_exptime(band,trange,
	#											verbose=verbose,retries=retries)

	#if len(tranges)==1:
	# Put the time range into the primary header for a single frame image
	#	hdu.header['EXPSTART'],hdu.header['EXPEND'] = tranges[0]
		# These are the proper keywords for this:
	#	hdu.header['TIME-OBS'],hdu.header['TIME-END'] = tranges[0]

	return hdu

def makemap(band,skypos,trange,skyrange,response=False,verbose=0):
	imsz = gxt.deg2pix(skypos,skyrange)
	photons = np.array(gQuery.getArray(gQuery.skyrect(band,
		skypos[0],skypos[1],trange[0],trange[1],skyrange[0],skyrange[1]),
		verbose=verbose),dtype='float64')
	try:
		events = {'t':photons[:,0 ]/tscale,'ra':photons[:,1],'dec':photons[:,2],
			'xi':photons[:,3],'eta':photons[:,4],
			'x':photons[:,5], 'y':photons[:,6]}
	except IndexError:
		if verbose>2:
			print 'No events found at {s} +/- {r} in {t}.'.format(
				s=skypos,r=skyrange,t=trange)
		return np.zeros(imsz)
	events = ct.hashresponse(band,events)
	wcs = define_wcs(skypos,skyrange,width=False,height=False)
	coo = zip(events['ra'],events['dec'])
	foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo,1),1)
	weights = 1./events['response'] if response else None
	H,xedges,yedges=np.histogram2d(foc[:,1]-0.5,foc[:,0]-0.5,bins=imsz,
		range=([ [0,imsz[0]],[0,imsz[1]] ]),weights=weights)
	return H

def integrate_map(band,skypos,tranges,skyrange,width=False,height=False,
				  verbose=0,memlight=False,hdu=False,retries=20,response=False):
	""" Integrate an image over some number of time ranges. Use a reduced
	memory optimization (at the expense of more web queries) if requested.
	"""
	imsz = gxt.deg2pix(skypos,skyrange)
	img = np.zeros(imsz)
	for trange in tranges:
		# If memlight is requested, break the integration into
		#  smaller chunks.
		# FIXME: memlight gives slightly wrong answers right now
		# This is probably due to a quirk of SQL, per issue #140.
		# Deprecating memlight until this can be resolved.
		step = memlight if memlight else trange[1]-trange[0]
		for i in np.arange(trange[0],trange[1],step):
			t0,t1=i,i+step
			if verbose:
				mc.print_inline('Coadding '+str(t0)+' to '+str(t1))
			img += makemap(band,skypos,[t0,t1],skyrange,response=response,
							 verbose=verbose)
		if response: # This is an intensity map.
			img /= dbt.compute_exptime(band,trange,skypos=skypos,
	                             						verbose=verbose)

	return img

def write_jpeg(filename,band,skypos,tranges,skyrange,width=False,height=False,
			   stepsz=1.,overwrite=False,verbose=0,retries=20):
	"""Write a 'preview' jpeg image from a count map."""
	scipy.misc.imsave(filename,integrate_map(band,skypos,tranges,skyrange,
					  width=width,height=height,verbose=verbose,
					  retries=retries))
	return

def movie(band,skypos,tranges,skyrange,framesz=0,width=False,height=False,
	verbose=0,memlight=False,coadd=False,response=False,hdu=False,retries=20):
	"""Generate a movie (mov)."""
	# Not defining stepsz creates a single full depth image.
	if coadd or (len(tranges)==1 and not framesz):
		if verbose>2:
			print 'Coadding across '+str(tranges)
		mv = integrate_map(band,skypos,tranges,skyrange,width=width,
			height=height,verbose=verbose,memlight=memlight,hdu=hdu,
			retries=retries,response=response)
		#rr.append(rrhr(band,skypos,tranges,skyrange,response=response,width=width,height=height,stepsz=1.,verbose=verbose,hdu=hdu,retries=retries)) if response else rr.append(np.ones(np.shape(mv)[1:]))
	else:
		for trange in tranges:
			stepsz = framesz if framesz else trange[1]-trange[0]
			steps = np.ceil((trange[1]-trange[0])/stepsz)
			for i,t0 in enumerate(np.arange(trange[0],trange[1],stepsz)):
				if verbose>1:
					mc.print_inline('Movie frame '+str(i+1)+' of '+
																str(int(steps)))
				t1 = trange[1] if i==steps else t0+stepsz

				img = integrate_map(band,skypos,[[t0,t1]],skyrange,
					width=width,height=height,verbose=verbose,
					memlight=memlight,hdu=hdu,retries=retries,
					response=response)
				if img.min() == 0 and img.max() == 0:
					if verbose>1:
						print 'No data in frame {i}. Skipping...'.format(i=i)
					continue
				try:
					mv.append(img)
				except:
					mv = [img]

	return np.array(mv)

def create_image(band,skypos,tranges,skyrange,framesz=0,width=False,
				 height=False,verbose=0,memlight=False,coadd=False,
				 response=False,hdu=False,retries=20):
	img = movie(band,skypos,tranges,skyrange,framesz=framesz,
		width=width,height=height,verbose=verbose,memlight=memlight,
		coadd=coadd,response=response,hdu=hdu,retries=retries)

	return np.array(img)

def write_images(band,skypos,tranges,skyrange,write_cnt=False,write_int=False,
				 write_rr=False,framesz=0,width=False,height=False,verbose=0,
				 memlight=False,coadd=False,overwrite=False,retries=20,
				 write_cnt_coadd=False, write_int_coadd=False):
	"""Generate a write various maps to files."""
	# No files were requested, so don't bother doing anything.
	imtypes = {'cnt':write_cnt,'int':write_int,'int_coadd':write_int_coadd,
			   'cnt_coadd':write_cnt_coadd}#,'rr':write_rr}
	for i in imtypes.keys():
		if not imtypes[i]:
			continue
		img = create_image(band,skypos,tranges,skyrange,framesz=framesz,
			width=width,height=height,verbose=verbose,memlight=memlight,
			retries=retries,
			coadd=True if (coadd or i in ['cnt_coadd','int_coadd']) else False,
			response=True if i in ['int','int_coadd'] else False)
		# Add a conditional so that this is only created for multi-frame images
		tbl = movie_tbl(band,tranges,framesz=framesz,verbose=verbose,
				retries=retries) if i in ['int','int_coadd'] else False
		hdu = pyfits.PrimaryHDU(img)
		hdu = fits_header(band,skypos,tranges,skyrange,width=width,
						  height=height,verbose=verbose,hdu=hdu,
						  retries=retries)
		hdulist = pyfits.HDUList([hdu,tbl]) if tbl else pyfits.HDUList([hdu])
		if verbose:
			print 'Writing image to {o}'.format(o=imtypes[i])
		hdulist.writeto(imtypes[i],clobber=overwrite)

	return
