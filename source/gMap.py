#!/usr/bin/python
import os
import ast
from imagetools import *
from dbasetools import fGetTimeRanges
from gphoton_args import setup_args, check_args

def gMap(band=None,calpath=os.pardir+os.sep+"cal"+os.sep,cntfile=False,coadd=False,detsize=1.25,intfile=False,rrfile=False,skypos=None,maxgap=1500.,memlight=100.,minexp=1.,overwrite=False,retries=20,skyrange=None,stepsz=0.,trange=None,verbose=0):
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

	# Use a mix of strings (if we want to make an output file) and Booleans (False if we do not).  I don't think this is the best way (I'd rather have separate variables for the True/False to create an image, and a string with the name of the output image), but for now this is kept since this is how the code was originally written.
	response = True if (intfile or rrfile) else False
	write_cnt = cntfile if (cntfile) else False
	write_int = intfile if (intfile) else False
	write_rr = rrfile if (rrfile) else False

	write_images(band,skypos,trange,skyrange,width=False,height=False,write_cnt=write_cnt,write_int=write_int,write_rr=write_rr,framesz=stepsz,clobber=overwrite,verbose=verbose,memlight=memlight,coadd=coadd,response=response,calpath=calpath,retries=retries)

if __name__ == '__main__':
	"""Called when gMap is executed directly through command line."""
	import argparse
	args = setup_args('gMap').parse_args()
	# Check command-line arguments.
	band,cntfile,detsize,intfile,rrfile,skypos,maxgap,memlight,minexp,retries,skyrange,stepsz,trange = check_args(args,"gMap")
	# Call gMap main program.
	gMap(band=band,calpath=args.calpath,cntfile=cntfile,coadd=args.coadd,detsize=detsize,intfile=intfile,rrfile=rrfile,skypos=skypos,maxgap=maxgap,memlight=memlight,minexp=minexp,overwrite=args.overwrite,retries=retries,skyrange=skyrange,stepsz=stepsz,trange=trange,verbose=args.verbose)
