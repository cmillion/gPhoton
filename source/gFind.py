#!/usr/bin/python
import ast
import argparse
import dbasetools as dbt
import gphoton_args as gargs

def gFind(band='both', detsize=1.25, exponly=False, gaper=False, maxgap=1.0,
		  minexp=1.0, quiet=False, retries=20, skypos=None, trange=None,
		  verbose=0,predicted=False):
	"""Primary program in the module. Prints time ranges to the screen and
    returns the total exposure time as a float.
	"""
	# Determine if we have to loop over both bands or just one.
	if band.upper() == 'BOTH':
		output = {'NUV':None,'FUV':None}
	elif band.upper() in ['NUV','FUV']:
		output = {band.upper():None}
	else:
		raise SystemExit('Invalid band: {b}'.format(b=band))
	all_expt = []
	for this_band in output.keys():
		## Get valid time ranges.
	       	ranges = dbt.fGetTimeRanges(this_band,skypos,maxgap=maxgap,
									minexp=minexp,trange=trange,verbose=verbose,
									detsize=detsize,retries=retries,
									predicted=predicted)

		if not len(ranges):
			if not quiet:
				print 'No {band} exposure time in database.'.format(
																band=this_band)
			output[this_band]={'expt':0,'t0':None,'t1':None}
		else:
			expt = (ranges[:,1]-ranges[:,0]).sum()
			if not quiet:
				print "{band}: {expt}s (raw) in {n} exposures.".format(
								band=this_band,expt=expt,n=len(ranges))
			if not exponly:
				if gaper:
					f = '['
					for r in ranges:
						f+='[%.3f' % r[0] + ', %.3f' % r[1] + '],'
					if not quiet: print f[:-1]+']'
				else:
					for r in ranges:
						if not quiet:
							print '    [ %.3f' % r[0] + ', %.3f' % r[1] + ' ], %.3f' % (r[1]-r[0]) + ' seconds'
			output[this_band]={'expt':expt,'t0':ranges[:,0],'t1':ranges[:,1]}
#			all_expt.append(expt)
	return output #{this_band:{'t0':ranges[:,0],'t1':ranges[:,1],'expt':all_expt}}

def setup_parser(iam='gfind'):
	parser = argparse.ArgumentParser(description="Locate available data.")
	parser = gargs.common_args(parser,iam)
	parser.add_argument("--alt", "--gaper", action="store_true",
		dest="gaper", default=False, help="Format the output so that it can "+
		"be copied and pasted directly into a gMap or gAperture command line?")
	parser.add_argument("--total", "--exponly", action="store_true",
		dest="exponly", default=False, help="Report only the total raw "+
		"exposure time available in the database.")
	parser.add_argument("--quiet", action="store_true", dest="quiet",
		help="Suppress all information to STDOUT.", default=False)
	parser.add_argument("--predicted", "--future", action="store_true",
		dest="predicted", help="Reported an estimate of the available"+
		"exposure _after the GR6/7 data is fully populated_.")
	return parser

def check_args(args,iam='gfind'):
	args = gargs.check_common_args(args,iam)
	return args

if __name__ == '__main__':
	"""Called when gFind is executed directly through command line."""
	args = setup_parser().parse_args()
	args = check_args(args)
	exp_times = gFind(band=args.band, detsize=args.detsize,
		exponly=args.exponly, gaper=args.gaper, maxgap=args.maxgap,
		minexp=args.minexp, quiet=args.quiet, retries=args.retries,
		skypos=args.skypos, trange=args.trange, predicted=args.predicted,
		verbose=args.verbose)
