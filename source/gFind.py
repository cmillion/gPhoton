#!/usr/bin/python

import ast
from dbasetools import fGetTimeRanges, suggest_parameters
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-b", "--band", action="store", type="string", dest="band", help="[NF]UV band designation", metavar="BAND")
parser.add_option("-r", "--ra", action="store", type="float", dest="ra", help="Center Right Ascension position in decimal degrees", metavar="RA")
parser.add_option("-d", "--dec", action="store", type="float", dest="dec", help="Center Declination position in decimal degrees", metavar="DEC")
parser.add_option("--skypos", action="store", type="string", dest="skypos", help="Alternate method for specifying sky position with format '[RA,Dec]'")
parser.add_option("-v", "--verbose", action="store", type="float", dest="verbose", help="Display more output. Set to 0-2", metavar="VRB",default=0)
parser.add_option("-g", "--gap", "--maxgap", action="store", type="float", dest="gap", help="Maximum gap size in seconds for data to be considered contiguous.", default=1)
parser.add_option("--t0", "--tmin", action="store", type="float", dest="tmin", help="Minimum time to consider.",default=1.)
parser.add_option("--t1", "--tmax", action="store", type="float", dest="tmax", help="Maxium time to consider.",default=1000000000000.)
parser.add_option("--trange", "--tranges", action="store", type="string", dest="trange", help="Time range in which to limit the search with format '[t0,t1]'")
parser.add_option("--detsize", action="store", type="float", dest="detsize", help="Set the effective field diameter in degrees for the exposure search.", default=1.25)
parser.add_option("--alt", "--gaper", action="store_true", dest="gaper", help="Format the output so that it can be copied and pasted directly into a gMap or gAperture command line.")
parser.add_option("--depth", "--exponly", action="store_true", dest="depth", help="Report only the total raw exposure time available in the database, not the individual time ranges.")
parser.add_option("--minexp", action="store", type="float", dest="minexp", help="Minimum contiguous exposure in seconds for data to be reported.", default=1)
parser.add_option("--suggest", action="store_true", dest="suggest", help="Suggest optimum parameters for aperture photometry.", default=False)

(options, args) = parser.parse_args()

mandatory = ['band']
for m in mandatory:
        if not options.__dict__[m]:
                print "A mandatory option is missing:",m
                parser.print_help()
                exit(1)

band = options.band
band = band.upper()
if band != "NUV" and band != "FUV":
        print "Band must be NUV or FUV. "
        exit(1)

if not (options.ra and options.dec) and not options.skypos:
	print "Must specify either RA/Dec or skypos."
	parser.print_help()
	exit(1)
elif (options.ra and options.dec) and options.skypos:
	print "Must specify either RA/Dec or skypos. Not both."
	parser.print_help()
	exit(1)
elif (options.ra and options.dec):
	skypos = [options.ra,options.dec]
else:
	skypos = list(ast.literal_eval(options.skypos))

if options.trange:
	trange = list(ast.literal_eval(options.trange))
else:
	trange = [options.tmin,options.tmax]
ranges = fGetTimeRanges(band,skypos,maxgap=options.gap,minexp=options.minexp,trange=trange,verbose=options.verbose,detsize=options.detsize)

if options.suggest:
	out = suggest_parameters(band,skypos,verbose=1)
	print ""

if not len(ranges):
	print 'No exposure time in database.'
else:
	expt   = (ranges[:,1]-ranges[:,0]).sum()
	print 'Available: '+str(expt)+' seconds (raw) in '+str(len(ranges))+' distinct exposures.'
	if not options.depth:
	        print ''
		if options.gaper:
			f = '['
			for r in ranges:
				f+='[%.3f' % r[0] + ', %.3f' % r[1] + '],'
			print f[:-1]+']'
		else:
			expt   = (ranges[:,1]-ranges[:,0]).sum()
			for r in ranges:
				print '    [ %.3f' % r[0] + ', %.3f' % r[1] + ' ], %.3f' % (r[1]-r[0]) + ' seconds'

