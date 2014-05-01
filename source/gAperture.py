#!/usr/bin/python

import ast
from curvetools import *
from imagetools import * # For JPEG preview image creation
from optparse import OptionParser
from dbasetools import fGetTimeRanges, suggest_parameters

parser = OptionParser()
parser.add_option("-b", "--band", action="store", type="string", dest="band", help="[NF]UV band designation")
parser.add_option("-r", "--ra", action="store", type="float", dest="ra", help="Center Right Ascension position in decimal degrees")
parser.add_option("-d", "--dec", action="store", type="float", dest="dec", help="Center Declination position in decimal degrees")
parser.add_option("--skypos", action="store", type="string", dest="skypos", help="Alternate method for specifying sky position with format '[RA,Dec]'")
parser.add_option("-a", "--aperture", action="store", type="float", dest="radius", help="Aperture radius in decimal degrees")
parser.add_option("-i", "--inner", action="store", type="float", dest="annulus1", help="Inner annulus radius for background subtraction")
parser.add_option("-o", "--outer", action="store", type="float", dest="annulus2", help="Outer annulus radius for background subtraction")
parser.add_option("--annulus", action="store", type="string", dest="annulus", help="Annulus inner and outer radius definition with format '[inner,outer]'")
parser.add_option("--t0", "--tmin", action="store", type="float", dest="tmin", help="Minimum time to consider.",default=1.)
parser.add_option("--t1", "--tmax", action="store", type="float", dest="tmax", help="Maxium time to consider.",default=1000000000000.)
parser.add_option("--tranges", action="store", type="string", dest="tranges", help="List of time ranges with format '[[t0,t1],[t2,t3],...]'")
parser.add_option("-s", "--step", "--frame", action="store", type="float", dest="stepsz", help="Step size in seconds",default=0)
parser.add_option("-f", "--file", action="store", type="string", dest="file", help="CSV output file")
parser.add_option("--response", action="store_true", dest="userr", help="Use the relative response correction.", default=False)
parser.add_option("--hrbg", action="store_true", dest="usehrbg", help="Use the higher quality 'swiss cheese' background estimation method.", default=False)
parser.add_option("-v", "--verbose", action="store", type="float", dest="verbose", help="Display more output. Set to 0-2.", default=0)
parser.add_option("--stamp", action="store", type="string", dest="stamp", help="Filename in which to write a preview JPEG image stamp of the targeted region.")
parser.add_option("--calpath", action="store", type="string", dest="calpath", help="Path to the directory that contains the calibration files.", default='../cal/')
parser.add_option("--addhdr", action="store_true", dest="addhdr", help="Write the command line and column headers to the top of the .csv file.")
parser.add_option("--coadd", action="store_true", dest="coadd", help="Return the coadded flux over all requested time ranges.")
parser.add_option("-g", "--gap", "--maxgap", action="store", type="float", dest="gap", help="Maximum gap size in seconds for data to be considered contiguous.", default=1)
parser.add_option("--minexp", action="store", type="float", dest="minexp", help="Minimum contiguous exposure in seconds for data to be reported.", default=1)
parser.add_option("--detsize", action="store", type="float", dest="detsize", help="Set the effective field diameter in degrees for the exposure search.", default=1.25)
parser.add_option("--bestparams", "--best", action="store_true", dest="best", help="Set parameters to produce the highest quality lightcurve. Potentially slow.", default=False)
parser.add_option("--suggest", "--optimize", action="store_true", dest="suggest", help="Suggest optimum parameters for aperture photometry.", default=False)

(options, args) = parser.parse_args()

# Reconstruct the command line
cmd = './gAperture.py'
for key in options.__dict__.keys():
	if options.__dict__[key]:
		cmd+=' --'+str(key)+' '+str(options.__dict__[key])

if not options.band:
	print "Band must be specified."
	exit(1)
else:
	options.band = options.band.upper()

if options.band != "NUV" and options.band != "FUV":
	print "Band must be NUV or FUV."
	exit(1)

if not options.radius and not options.suggest:
	print "Must specify an aperture radius."
	exit(1)

if not (options.ra and options.dec) and not options.skypos:
	print "Must specify either RA/Dec or skypos."
	exit(1)
elif (options.ra and options.dec) and options.skypos:
	print "Must specify either RA/Dec or skypos. Not both."
	exit(1)
elif (options.ra and options.dec):
	skypos = [options.ra,options.dec]
else:
	skypos = list(ast.literal_eval(options.skypos))

if options.suggest:
	options.ra,options.dec,options.radius,options.annulus1,options.annulus2 = suggest_parameters(options.band,skypos)
	if options.verbose:
		print "Recentering on ["+str(options.ra)+", "+str(options.dec)+"]"
		print "Setting radius to "+str(options.radius)
		print "Setting annulus to ["+str(options.annulus1)+", "+str(options.annulus2)+"]"

if options.stepsz and options.coadd:
	print "Cannot specify both --stepsz and --coadd."
	exit(1)

if options.best:
	options.userr = True
	options.usehrbg = True

# This breaks the 'high resolution' background
if not (options.annulus1 and options.annulus2) and not options.annulus:
	if options.usehrbg:
		print "Must specify background annulus with hrbg."
		exit(1)
	else:
		annulus = [False,False]
elif options.annulus1 and options.annulus2:
	annulus = [options.annulus1,options.annulus2]
else:
	annulus = ast.literal_eval(options.annulus)

if options.tranges:
	tranges = list(ast.literal_eval(options.tranges))
else:
	if options.verbose:
		print 'Using all available exposure time between '+str(options.tmin)+' and '+str(options.tmax)
	tranges = fGetTimeRanges(options.band,skypos,maxgap=options.gap,verbose=options.verbose,minexp=options.minexp,trange=[options.tmin,options.tmax],detsize=options.detsize)
	if not len(tranges):
		print 'No exposure time in database.'
		exit(0)
	if options.verbose:
		print 'Using '+str((tranges[:,1]-tranges[:,0]).sum())+' seconds (raw) in '+str(len(tranges))+' distinct exposures.'

# Make sure tranges is a 2D array
if len(np.array(tranges).shape)==1:
	tranges=[tranges]

iocode = 'wb'
if not options.file:
	outfile = False
else:
	# This initiates the file.
	outfile = options.file
	f = open(outfile,iocode)
	if options.addhdr:
		# Write the command line to the outfile; should be optional
		f.write('| '+cmd+'\n')
		f.write('| tstart, tstop, radius, exptime, cps, error, flux, flux_error, magnitude, mag_error, inner annulus, outer annulus, background, response, counts, aperture correction 1, aperture correction 2\n')
	f.close()
	# Setting this will now append to the file we just created
	iocode = 'ab'

if options.stamp:
	if options.annulus2:
		skyrange = [2*options.annulus2,2*options.annulus2]
	else:
		skyrange = [2*options.radius,2*options.radius]
	
	write_jpeg(options.stamp,options.band,skypos,tranges,skyrange,width=False,height=False,stepsz=options.stepsz,clobber=True,verbose=options.verbose)

write_curve(options.band,skypos,tranges,options.radius,outfile=outfile,annulus=annulus,stepsz=options.stepsz,userr=options.userr,usehrbg=options.usehrbg,verbose=options.verbose,calpath=options.calpath,iocode=iocode,coadd=options.coadd,detsize=options.detsize)
