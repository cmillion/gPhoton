import csv
from astropy.io import fits as pyfits
import numpy as np
from sys import stdout

import pyximport
pyximport.install()

# A package for html queries
import requests
from MCUtils import *
import gQuery

def load_raw6(raw6file):
	"""Reads a raw6 file. Just wraps some pyfits commands."""
	print "		",raw6file
	hdulist = pyfits.open(raw6file,memmap=1)
	htab = hdulist[1].header
	hdulist.close()
	return htab, hdulist

def find_band(raw6file):
	"""Tries to guess the band based on a raw6 filename."""
	# This assumes the standard raw6 naming convention. I made it
	# this needlessly specific to avoid false positives in the directory path.
	if raw6file.rfind('-fd-raw6.fits') > 0 and raw6file.rfind('-nd-raw6.fits') > 0:
		print "Multiple bands implied in raw6 filename. Specify band on command line."
		return 0

	if raw6file.rfind('-fd-raw6.fits') > 0:
		band = 'FUV'
	elif raw6file.rfind('-nd-raw6.fits') > 0:
		band = 'NUV'
	else:
		print "Band not derivable from raw6 filename. Specify band on command line."
		band = 0

	return band

def load_aspect(aspfile):
	"""Loads an aspect file into a bunch of arrays."""
	ra,dec,twist,time,aspflags=np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
	header = {'RA':[],'DEC':[],'ROLL':[]}
	for i in xrange(len(aspfile)):
        	print "         ",aspfile[i]
		hdulist = pyfits.open(aspfile[i],memmap=1)
		ra = np.append(ra,np.array(hdulist[1].data.field('ra')))
		dec = np.append(dec,np.array(hdulist[1].data.field('dec')))
		twist = np.append(twist,np.array(hdulist[1].data.field('roll')))
		time = np.append(time,np.array(hdulist[1].data.field('t')))
		aspflags = np.append(aspflags,np.array(hdulist[1].data.field('status_flag')))
		header['RA'] = np.append(header['RA'],np.zeros(len(hdulist[1].data.field('ra')))+hdulist[0].header['RA_CENT'])
		header['DEC'] = np.append(header['DEC'],np.zeros(len(hdulist[1].data.field('dec')))+hdulist[0].header['DEC_CENT'])
		header['ROLL'] = np.append(header['ROLL'],np.zeros(len(hdulist[1].data.field('roll')))+hdulist[0].header['ROLL'])
	#header = {'RA':np.array(ra).mean(),'DEC':np.array(dec).mean(),'ROLL':np.array(twist).mean()}

	# Time sort stuff.
	ix = np.argsort(time)

	hdulist.close()
	#header['RA'] = header['RA'][ix]
	#header['DEC'] = header['DEC'][ix]
	#header['ROLL'] = header['ROLL'][ix]

	return ra[ix], dec[ix], twist[ix], time[ix], header, aspflags[ix]

def web_query_aspect(eclipse):
	"""Grabs the aspect data from MAST databases based on eclipse."""
	print "Attempting to query MAST database for aspect records."
	entries = gQuery.getArray(gQuery.aspect_ecl(eclipse))
	n = len(entries)
	print '		Located '+str(n)+' aspect entries.'
	if not n:
		print "No aspect entries for eclipse "+str(eclipse)
		return
	ra,dec,twist,time,flags=[],[],[],[],[]
	header = {'RA':[],'DEC':[],'ROLL':[]}
	ra0,dec0,twist0 = [],[],[]
	for i in xrange(n):
		# The times are *1000 in the database to integerify
		time.append(float(entries[i][2])/1000.)
		ra.append(float(entries[i][3]))
		dec.append(float(entries[i][4]))
		twist.append(float(entries[i][5]))
		flags.append(float(entries[i][6]))
		ra0.append(float(entries[i][7]))
		dec0.append(float(entries[i][8]))
		twist0.append(float(entries[i][9]))

	# Need to sort the output so that it is time ordered before returning.
	# Although it should already be ordered by time because that is requested
	#  in the SQL query above. If this is time consuming, remove it.
	ix = np.argsort(np.array(time))
	header = {'RA':np.array(ra0)[ix],'DEC':np.array(dec0)[ix],'ROLL':np.array(twist0)[ix]}

	return np.array(ra)[ix],np.array(dec)[ix],np.array(twist)[ix],np.array(time)[ix],header,np.array(flags)[ix]

def wiggle_filenames(band,calpath):
	"""Returns the 'wiggle' calibration file name."""
	if band == 'NUV':
		wiggle_files = {'x':calpath+'NUV_wiggle_x.fits','y':calpath+'NUV_wiggle_y.fits'}
	elif band == 'FUV':
		wiggle_files = {'x':calpath+'FUV_wiggle_x.fits','y':calpath+'FUV_wiggle_y.fits'}
	else:
		print "Band not specified."

	return wiggle_files

def avgwalk_filenames(band,calpath):
	"""Returns the 'avgwalk' calibration filename."""
	if band == 'NUV':
		walk_files = {'x':calpath+'NUV_avgwalk_x.fits','y':calpath+'NUV_avgwalk_y.fits'}
	elif band == 'FUV':
		walk_files = {'x':calpath+'FUV_avgwalk_x.fits','y':calpath+'FUV_avgwalk_y.fits'}
	else:
		print "Band not specified."

	return walk_files

def walk_filenames(band,calpath):
	"""Returns the 'walk' calibration filename."""
	if band == 'NUV':
		walk_files = {'x':calpath+'NUV_walk_x.fits','y':calpath+'NUV_walk_y.fits'}
	elif band == 'FUV':
		walk_files = {'x':calpath+'FUV_walk_x.fits','y':calpath+'FUV_walk_y.fits'}
	else:
		print "Band not specified."

	return walk_files

def linearity_filenames(band,calpath):
	"""Returns the 'linearity' calibration file name."""
	if band == 'NUV':
		linfiles = {'x':calpath+'NUV_NLC_x_det2sky.fits','y':calpath+'NUV_NLC_y_det2sky.fits'}
	elif band == 'FUV':
		linfiles = {'x':calpath+'FUV_NLC_x_det2sky.fits','y':calpath+'FUV_NLC_y_det2sky.fits'}
	else:
		print "Band not specified."

	return linfiles

def flat_filename(band,calpath):
	"""Returns the 'flat' calibration file name."""
	if band=='NUV':
		flatfile = calpath+'NUV_flat.fits'
	elif band=='FUV':
		flatfile = calpath+'FUV_flat.fits'
	else:
		print "Band not specified."

	return flatfile

def distortion_filenames(band,calpath,eclipse,raw_stimsep):
	"""Returns the 'distortion' calibration file names."""
	if band == 'NUV':
		if (eclipse > 37460):
			if (raw_stimsep < 5136.3):
				distfiles = {'x':calpath+'nuv_distortion_cube_dxa.fits','y':calpath+'nuv_distortion_cube_dya.fits'}
			elif (raw_stimsep < 5137.25):
				distfiles = {'x':calpath+'nuv_distortion_cube_dxb.fits','y':calpath+'nuv_distortion_cube_dyb.fits'}
			else:
				distfiles = {'x':calpath+'nuv_distortion_cube_dxc.fits','y':calpath+'nuv_distortion_cube_dyc.fits'}
		else:
			distfiles = {'x':calpath+'nuv_distortion_cube_dx.fits','y':calpath+'nuv_distortion_cube_dy.fits'}
	elif band == 'FUV':
		distfiles = {'x':calpath+'fuv_distortion_cube_dx.fits','y':calpath+'fuv_distortion_cube_dy.fits'}
	else:
		print "Band not specified."

	return distfiles

def offset_filenames(calpath):
	"""Returns the NUV->FUV offset calibration file names."""
	# Offset is for FUV only
	return {'x':calpath+'fuv_dx_fdttdc_coef_0.tbl','y':calpath+'fuv_dy_fdttdc_coef_0.tbl'}

def mask_filename(band,calpath):
	if band == 'NUV':
		maskfile = calpath+'NUV_mask.fits'
	elif band == 'FUV':
		maskfile = calpath+'FUV_mask.fits'
	else:
		print "Band not specified."

	return maskfile

def create_SSD_filename(band,eclipse):
	"""Returns the Stim Separation Data (SSD) calibration file name."""
	return "SSD_"+band.lower()+"_"+str(eclipse)+".tbl"

