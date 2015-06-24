import csv
from astropy.io import fits as pyfits
import numpy as np
from sys import stdout

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

def web_query_aspect(eclipse,retries=20):
	"""Grabs the aspect data from MAST databases based on eclipse."""
	print "Attempting to query MAST database for aspect records."
	entries = gQuery.getArray(gQuery.aspect_ecl(eclipse),retries=retries)
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

def create_SSD_filename(band,eclipse):
	"""Returns the Stim Separation Data (SSD) calibration file name."""
	return "SSD_"+band.lower()+"_"+str(eclipse)+".tbl"
