#!/usr/bin/python
import os
import ast
import sys
from curvetools import *
from imagetools import * # For JPEG preview image creation
from optparse import OptionParser
from dbasetools import fGetTimeRanges, suggest_parameters
from gphoton_args import setup_args, check_args

def gAperture(addhdr=False,annulus=None,annulus1=None,annulus2=None,band=None,best=False,calpath=os.pardir+os.sep+"cal"+os.sep,cmd=None,coadd=False,dec=None,detsize=1.25,outfile=None,maxgap=1500.,minexp=1.,overwrite=False,ra=None,radius=None,retries=20,skypos=None,stamp=None,stepsz=0.,suggest=False,t0=1.,t1=1000000000000.,trange=None,usehrbg=False,userr=False,verbose=0):
	"""Primary program in the module.  Generates lightcurve CSV files and preview stamp images."""

	if verbose:
		print 'Using all exposure in ['+str(trange[0])+','+str(trange[1])+']'
	trange = fGetTimeRanges(band,skypos,maxgap=maxgap,verbose=verbose,minexp=minexp,trange=trange,detsize=detsize,retries=retries)
	if not len(trange):
		print 'No exposure time in database.'
		exit(0)
	if verbose:
		print 'Using '+str((trange[:,1]-trange[:,0]).sum())+' seconds (raw) in '+str(len(trange))+' distinct exposures.'

	# Make sure trange is a 2D array
	if len(np.array(trange).shape)==1:
		trange=[trange]

	iocode = 'wb'
	if not outfile:
		outfile = False
	else:
		# This initiates the file.
		if not overwrite and os.path.exists(outfile):
			print "File "+str(outfile)+"exists.  Use --overwrite to replace it."
			exit(0)
		f = open(outfile,iocode)
		if addhdr:
			# Write the command line to the outfile; should be optional
			f.write('| '+cmd+'\n')
			f.write('| tstart, tstop, radius, exptime, cps, error, flux, flux_error, magnitude, mag_error, inner annulus, outer annulus, background, response, counts, aperture correction 1, aperture correction 2\n')
			f.close()
			# Setting this will now append to the file we just created
			iocode = 'ab'

	if stamp:
		if annulus:
			skyrange = [2.*annulus[1],2.*annulus[1]]
		else:
			skyrange = [2.*radius,2.*radius]
		write_jpeg(stamp,band,skypos,trange,skyrange,width=False,height=False,stepsz=stepsz,clobber=overwrite,verbose=verbose,retries=retries)

	write_curve(band,skypos,trange,radius,outfile=outfile,annulus=annulus,stepsz=stepsz,userr=userr,usehrbg=usehrbg,verbose=verbose,calpath=calpath,iocode=iocode,coadd=coadd,detsize=detsize,retries=retries)

if __name__ == '__main__':
	"""Called when gAperture is executed directly through command line."""
	import argparse
	args = setup_args('gAperture').parse_args()
	# Reconstruct the command line to pass to gAperture if ADDHDR is requested.
	cmd = " ".join(sys.argv)
	# Check command-line arguments.
	annulus,band,best,detsize,outfile,maxgap,minexp,radius,retries,skypos,stamp,stepsz,trange,usehrbg,userr = check_args(args,"gAperture")
	# Call gAperture main program.
	gAperture(addhdr=args.addhdr,annulus=annulus,band=band,best=best,calpath=args.calpath,cmd=cmd,coadd=args.coadd,detsize=detsize,outfile=outfile,maxgap=maxgap,minexp=minexp,overwrite=args.overwrite,radius=radius,retries=retries,skypos=skypos,stamp=stamp,stepsz=stepsz,suggest=args.suggest,trange=trange,usehrbg=usehrbg,userr=userr,verbose=args.verbose)
