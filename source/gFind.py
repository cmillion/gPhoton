#!/usr/bin/python

import ast
from dbasetools import fGetTimeRanges, suggest_parameters
from gphoton_args import setup_args, check_args

def gFind(band='both', detsize=1.25, exponly=False, gaper=False, maxgap=1.0, minexp=1.0, quiet=False, retries=20, skypos=None, trange=None, verbose=0):
	"""Primary program in the module. Prints time ranges to the screen and returns
	the total exposure time as a float.
	"""
	# Determine if we have to loop over both bands or just one.
	if band.lower() == 'both':
		all_bands = ['fuv','nuv']
	else:
		all_bands = [band]
	all_expt = []

	for this_band in all_bands:
		## Get valid time ranges.
	       	ranges = fGetTimeRanges(this_band,skypos,maxgap=maxgap,minexp=minexp,trange=trange,verbose=verbose,detsize=detsize,retries=retries)

		if not len(ranges):
			if not quiet: print 'No '+this_band.upper()+' exposure time in database.'
			return 0.
		else:
			expt = (ranges[:,1]-ranges[:,0]).sum()
			if not quiet: print this_band.upper()+': Available: '+str(expt)+' seconds (raw) in '+str(len(ranges))+' distinct exposures.'
			if not exponly:
				if not quiet: print ''
				if gaper:
					f = '['
					for r in ranges:
						f+='[%.3f' % r[0] + ', %.3f' % r[1] + '],'
					if not quiet: print f[:-1]+']'
				else:
					for r in ranges:
						if not quiet: print '    [ %.3f' % r[0] + ', %.3f' % r[1] + ' ], %.3f' % (r[1]-r[0]) + ' seconds'
			all_expt.append(expt)
	return all_expt

if __name__ == '__main__':
	"""Called when gFind is executed directly through command line."""
	import argparse
	args = setup_args('gFind').parse_args()
	band, detsize, maxgap, minexp, retries, skypos, trange = check_args(args,"gFind")
	exp_times = gFind(band=band, detsize=detsize, exponly=args.exponly, gaper=args.gaper, maxgap=maxgap, minexp=minexp, quiet=args.quiet, retries=retries, skypos=skypos, trange=trange, verbose=args.verbose)
