#!/usr/bin/python

import ast
from dbasetools import fGetTimeRanges, suggest_parameters

def gFind(band, skypos, trange, maxgap=1.0, minexp=1.0, verbose=0., detsize=1.25, depth=False, gaper=False, quiet=False):
	"""Primary program module. Prints time ranges to the screen and returns
	the total exposure time as a float.
	"""
	## Get valid time ranges.
       	ranges = fGetTimeRanges(band,skypos,maxgap=maxgap,minexp=minexp,trange=trange,verbose=verbose,detsize=detsize)

	if not len(ranges):
		if not quiet: print 'No exposure time in database.'
		return 0.
	else:
		expt = (ranges[:,1]-ranges[:,0]).sum()
		if not quiet: print 'Available: '+str(expt)+' seconds (raw) in '+str(len(ranges))+' distinct exposures.'
		if not depth:
			if not quiet: print ''
			if gaper:
				f = '['
				for r in ranges:
					f+='[%.3f' % r[0] + ', %.3f' % r[1] + '],'
				if not quiet: print f[:-1]+']'
			else:
				for r in ranges:
					if not quiet: print '    [ %.3f' % r[0] + ', %.3f' % r[1] + ' ], %.3f' % (r[1]-r[0]) + ' seconds'
		return expt

def check_args(args):
	"""Checks validity of some command line arguments."""
	# Make sure "band" is defined and is a valid choice.
	mandatory = ['band']
	for m in mandatory:
		if not args.__dict__[m]:
			print "A mandatory option is missing:",m
			exit(1)
	band = args.band
	band = band.upper()
	if band != "NUV" and band != "FUV":
		print "Band must be NUV or FUV. "
		exit(1)

	# Create the "skypos" list whether or not skypos was directly submitted
	# or the RA and DEC were given separately.
	if not (args.ra and args.dec) and not args.skypos:
		print "Must specify either RA/Dec or skypos."
		exit(1)
	elif (args.ra and args.dec) and args.skypos:
		print "Must specify either RA/Dec or skypos. Not both."
		exit(1)
	elif (args.ra and args.dec):
		skypos = [args.ra,args.dec]
	else:
		skypos = list(ast.literal_eval(args.skypos))

	# Print suggested parameters to screen if requested.
	if args.suggest:
		out = suggest_parameters(band,skypos,verbose=1)
		print ""

	# Create the "trange" list whether it's given directly or
	# if tmin and tmax are specified separately.
       	if args.trange:
		trange = list(ast.literal_eval(args.trange))
	else:
		trange = [args.tmin,args.tmax]

	# Make sure the "gap" parameter is greater or equal to 1.
	if args.gap >= 1.0:
		maxgap = args.gap
	else:
		print "Maximum gap size must be >= 1.0."
		exit(1)

	# Make sure the "minexp" parameter is greater or equal to 1.
	if args.minexp >= 1.0:
		minexp = args.minexp
	else:
		print "Minimum contiguous exposure length must be >= 1.0."
		exit(1)

	# Make sure the "detsize" parameter is greater than 0.
	if args.detsize > 0.0:
		detsize = args.detsize
	else:
		print "Minimum effective field diameter must be > 0.0."
		exit(1)

	return (band, skypos, trange, maxgap, minexp, detsize)

def setup_parser():
	"""Defines the arguments and options for the parser object when
	called from the command line.
	"""
	parser = argparse.ArgumentParser(description="Locate available data at specified coordinates and time intervals.")
	parser.add_argument("-b", "--band", action="store", dest="band", help="[NF]UV band designation", metavar="BAND")
	parser.add_argument("-r", "--ra", action="store", type=float, dest="ra", help="Center Right Ascension position in decimal degrees", metavar="RA")
	parser.add_argument("-d", "--dec", action="store", type=float, dest="dec", help="Center Declination position in decimal degrees", metavar="DEC")
	parser.add_argument("--skypos", action="store", dest="skypos", help="Alternate method for specifying sky position with format '[RA,Dec]'")
	parser.add_argument("-v", "--verbose", action="store", type=float, dest="verbose", help="Display more output. Set to 0-2", metavar="VRB",default=0.,choices=[0,1,2])
	parser.add_argument("-g", "--gap", "--maxgap", action="store", type=float, dest="gap", help="Maximum gap size in seconds for data to be considered contiguous.", default=1.)
	parser.add_argument("--t0", "--tmin", action="store", type=float, dest="tmin", help="Minimum time to consider.",default=1.)
	parser.add_argument("--t1", "--tmax", action="store", type=float, dest="tmax", help="Maxium time to consider.",default=1000000000000.)
	parser.add_argument("--trange", "--tranges", action="store", dest="trange", help="Time range in which to limit the search with format '[t0,t1]'")
	parser.add_argument("--detsize", action="store", type=float, dest="detsize", help="Set the effective field diameter in degrees for the exposure search.", default=1.25)
	parser.add_argument("--alt", "--gaper", action="store_true", dest="gaper", help="Format the output so that it can be copied and pasted directly into a gMap or gAperture command line.", default=False)
	parser.add_argument("--depth", "--exponly", action="store_true", dest="depth", help="Report only the total raw exposure time available in the database, not the individual time ranges.", default=False)
	parser.add_argument("--minexp", action="store", type=float, dest="minexp", help="Minimum contiguous exposure in seconds for data to be reported.", default=1.)
	parser.add_argument("--suggest", action="store_true", dest="suggest", help="Suggest optimum parameters for aperture photometry.", default=False)
	parser.add_argument("--quiet", action="store_true", dest="quiet", help="Suppress all information to STDOUT.", default=False)
	return parser

if __name__ == '__main__':
	"""Called when gFind is executed directly through command line."""
	import argparse
	args = setup_parser().parse_args()
	band, skypos, trange, maxgap, minexp, detsize = check_args(args)
	gFind(band, skypos, trange, maxgap, minexp, args.verbose, detsize, args.depth, args.gaper,args.quiet)

