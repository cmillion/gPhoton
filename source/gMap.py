#!/usr/bin/python
import os
import ast
from imagetools import *
from dbasetools import fGetTimeRanges
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--count", action="store", type="string", dest="cntfile", help="Filename for the count image.", default=False)
parser.add_option("--intensity", action="store", type="string", dest="intfile", help="Filename for the intensity image.", default=False)
parser.add_option("--response", action="store", type="string", dest="rrfile", help="Filename for the response image.", default=False)
parser.add_option("-b", "--band", action="store", type="string", dest="band", help="[NF]UV band designation")
parser.add_option("-r", "--ra", action="store", type="float", dest="ra", help="Center Right Ascension position in decimal degrees")
parser.add_option("-d", "--dec", action="store", type="float", dest="dec", help="Center Declination position in decimal degrees")
parser.add_option("--skypos", action="store", type="string", dest="skypos", help="Alternate method for specifying sky position with format '[RA,Dec]'")
parser.add_option("--angle", action="store", type="float", dest="angle", help="The angle subtended in both RA and Dec.")
parser.add_option("--raangle", action="store", type="float", dest="raangle", help="The angle of sky in degrees that the right ascension subtends. Overrides --angle")
parser.add_option("--decangle", action="store", type="float", dest="decangle", help="The angle of sky in degrees that the declination subtends. Overrides --angle")
parser.add_option("--t0", "--tmin", action="store", type="float", dest="tmin", help="Minimum time to consider.",default=1.)
parser.add_option("--t1", "--tmax", action="store", type="float", dest="tmax", help="Maxium time to consider.",default=1000000000000.)
parser.add_option("--trange", action="store", type="string", dest="trange", help="List of time ranges with format '[[t0,t1,],[t2,t3],...]'")
parser.add_option("--frame", "--step", action="store", type="float", dest="stepsz", help="Depth of movie frames in seconds.", default=0.)
parser.add_option("--calpath", action="store", type="string", dest="calpath", help="Path to the directory that contains the calibration files.", default='../cal/')
parser.add_option("-v", "--verbose", action="store", type="float", dest="verbose", help="Display more output. Set to 0-2.", default=0)
parser.add_option("--memlight", action="store", type="float", dest="memlight", help="Reduce server-side memory usage by requesting data in chunks of no more than this depth in seconds.", default=100)
parser.add_option("--coadd", action="store_true", dest="coadd", help="Return an image coadded over all requested time ranges.")
parser.add_option("-g", "--gap", "--maxgap", action="store", type="float", dest="gap", help="Maximum gap size in seconds for data to be considered contiguous.", default=1)
parser.add_option("--minexp", action="store", type="float", dest="minexp", help="Minimum contiguous exposure in seconds for data to be reported.", default=1)
parser.add_option("--maskrad", action="store", type="float", dest="maskrad", help="The radius at which detector events will be masked out.", default=1.)
parser.add_option("--overwrite", "--ow", "--clobber", action="store_true", dest="overwrite", default=False)
parser.add_option("--retries", action="store", type="int", dest="retries", default=20, help="Set the number of times to ping the server for a response before defining a query failure.  Default is 20, set to a large number if you expect, or want to allow, the query to take a long time.")

(options, args) = parser.parse_args()

mandatory = ['band']
for m in mandatory:
	if not options.__dict__[m]:  
		print "A mandatory option is missing:",m
		parser.print_help()
		exit(1)

options.band = options.band.upper()
if options.band != "NUV" and options.band != "FUV":
	print "Band must be NUV or FUV."
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

if options.stepsz and options.coadd:
	print "Cannot specify by --frame and --coadd."
	exit(1)

# Default to using memlight if this is NUV
if not options.memlight and options.band=="NUV":
	options.memlight = 100.

if options.angle:
	skyrange = [options.angle,options.angle]
	if options.raangle:
		skyrange[0] = options.raangle
	if options.decangle:
		skyrange[0] = options.decangle
elif options.raangle and options.decangle:
	skyrange = [options.raangle,options.decangle]
else:
	print "Must specify either --raangle and --decangle or --angle."
	exit(1)

if not options.overwrite:
	if options.cntfile and os.path.exists(options.cntfile):
		print 'File '+str(options.cntfile)+' exists. Use --overwrite.'
		exit(0)
	if options.cntfile and os.path.exists(options.cntfile):
		print 'File '+str(options.cntfile)+' exists. Use --overwrite.'
		exit(0)
	if options.cntfile and os.path.exists(options.cntfile):
		print 'File '+str(options.cntfile)+' exists. Use --overwrite.'
		exit(0)

if options.trange:
	trange = list(ast.literal_eval(options.trange))
else:
	if options.verbose and not (options.tmin or options.tmax):
		print 'No time range given. Using all available exposure time.'
	else:
		print 'Using all exposure in ['+str(options.tmin)+','+str(options.tmax)+']'
	trange = fGetTimeRanges(options.band,skypos,maxgap=options.gap,verbose=options.verbose,minexp=options.minexp,trange=[options.tmin,options.tmax],retries=options.retries)
	if not len(trange):
		print 'No exposure time in database.'
		exit(0)
	if options.verbose:
		print 'Using '+str((trange[:,1]-trange[:,0]).sum())+' seconds (raw) in '+str(len(trange))+' distinct exposures.'

# Make sure trange is a 2D array
if len(np.array(trange).shape)==1:
	trange=[trange]


response = True if (options.intfile or options.rrfile) else False

write_images(options.band,skypos,trange,skyrange,width=False,height=False,write_cnt=options.cntfile,write_int=options.intfile,write_rr=options.rrfile,framesz=options.stepsz,clobber=options.overwrite,verbose=options.verbose,memlight=options.memlight,coadd=options.coadd,response=response,calpath=options.calpath,retries=options.retries)

