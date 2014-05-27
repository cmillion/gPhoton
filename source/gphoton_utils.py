# read and plot functionality for gPhoton .csv lightcurve files as created by gAperture.
import os
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from galextools import counts2mag, mag2counts

def read_lc(ifile,header=False):
    """Reads a light curve CSV into a pandas DataFrame"""
    ## Read in the data as a pandas DataFrame object.
    try:
        data_frame = pd.io.parsers.read_csv(ifile, skipinitialspace=True, names=["t0", "t1", "ap_radius", "exptime", "cps", "cpserr", "flux", "flux_err", "mag", "mag_err", "r_inner", "r_outer", "bkg", "response", "counts", "apcorr1", "apcorr2"],skiprows=(2 if header else False))
        ## Calculate the timestamp in JD.
        data_frame["JD"] = ((data_frame["t1"] - data_frame["t0"]) / 2. + data_frame["t0"] + 315964800.0) / 86400. + 2440587.5
        
        ## Replace the default index with the timestamp in JD.
        data_frame.set_index("JD", drop=True, inplace=True, verify_integrity=True)

        ## Return the pandas objects.
        return data_frame
    except IOError:
        print "*** Error:  Could not find the file " + ifile + "."
	return

def plot_lc(data_frame):
	"""Plots a lightcurve from a CSV file
	data_frame - pandas DataFrame from read_lc()
	"""
	plt.plot(data_frame.index.values, data_frame["flux"], "ko")
	plt.show()
	return

def model_errors(catmag,band,sigma=3,mode='mag',trange=[1,1600]):
	"""Give upper and lower expected bounds as a function of the nominal
	magnitude of a source. Super userful for identifying outliers.
	sigma = how many sigma out to set the bounds
	mode = ['cps','mag'] - units in which to report bounds
	"""
	if mode!='cps' and mode!='mag':
		print 'mode must be set to "cps" or "mag"'
		exit(0)
	x = np.arange(trange[0],trange[1])
	cnt = mag2counts(catmag,band)
	ymin = (cnt*x/x)-sigma*np.sqrt(cnt*x)/x
	ymax = (cnt*x/x)+sigma*np.sqrt(cnt*x)/x
	if mode=='mag':
#		y = counts2mag(cnt*x/x,band)
		ymin = counts2mag(ymin,band)
		ymax = counts2mag(ymax,band)
	return ymin, ymax

# Given an array (of counts or mags) return an array of 1-sigma error values
def data_errors(catmag,t,band,sigma=3,mode='mag'):
	if mode!='cps' and mode!='mag':
		print 'mode must be set to "cps" or "mag"'
		exit(0)
	cnt = mag2counts(catmag,band)
	ymin = (cnt*t/t)-sigma*np.sqrt(cnt*t)/t
	ymax = (cnt*t/t)+sigma*np.sqrt(cnt*t)/t
	if mode=='mag':
#		y = counts2mag(cnt*t/t,band)
		ymin = counts2mag(ymin,band)
		ymax = counts2mag(ymax,band)
	return ymin, ymax


