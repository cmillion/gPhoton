# read and plot functionality for gPhoton .csv lightcurve files as created by gAperture.
import os
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import galextools as gt

def read_lc(csvfile,comment='|'):
	return pd.io.parsers.read_csv(csvfile,comment=comment)

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
	cnt = gt.mag2counts(catmag,band)
	ymin = (cnt*x/x)-sigma*np.sqrt(cnt*x)/x
	ymax = (cnt*x/x)+sigma*np.sqrt(cnt*x)/x
	if mode=='mag':
#		y = gt.counts2mag(cnt*x/x,band)
		ymin = gt.counts2mag(ymin,band)
		ymax = gt.counts2mag(ymax,band)
	return ymin, ymax

# Given an array (of counts or mags) return an array of 1-sigma error values
def data_errors(catmag,t,band,sigma=3,mode='mag'):
	if mode!='cps' and mode!='mag':
		print 'mode must be set to "cps" or "mag"'
		exit(0)
	cnt = gt.mag2counts(catmag,band)
	ymin = (cnt*t/t)-sigma*np.sqrt(cnt*t)/t
	ymax = (cnt*t/t)+sigma*np.sqrt(cnt*t)/t
	if mode=='mag':
#		y = gt.counts2mag(cnt*t/t,band)
		ymin = gt.counts2mag(ymin,band)
		ymax = gt.counts2mag(ymax,band)
	return ymin, ymax

def dmag_errors(t,band,sigma=3,mode='mag',mags=np.arange(13,24,0.1)):
    """Given an exposure time, give dmag error bars at a range of magnitudes."""
    cnts = gt.mag2counts(mags,band)
    ymin = (cnts*t/t)-sigma*np.sqrt(cnts*t)/t
    ymax = (cnts*t/t)+sigma*np.sqrt(cnts*t)/t
    if mode=='mag':
        ymin=mags-gt.counts2mag(ymin,band)
        ymax=mags-gt.counts2mag(ymax,band)
    return mags,ymin,ymax
