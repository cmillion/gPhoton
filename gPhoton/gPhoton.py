#!/usr/bin/python

import os

from PhotonPipe import *

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--raw6", action="store", type="string", dest="raw6file", help="raw6 path/filename", metavar="FILE")
parser.add_option("-a", "--aspect", action="store", type="string", dest="aspfile", help="aspection path/filename", metavar="FILE")
parser.add_option("-s", "--scst", action="store", type="string", dest="scstfile", help="spacecraft state path/filename", metavar="SCST")
parser.add_option("-c", "--calpath", action="store", type="string", dest="calpath", help="calibration directory path", metavar="PATH")
parser.add_option("-b", "--band", action="store", type="string", dest="band", help="[NF]UV band designation", metavar="BAND")
parser.add_option("-o", "--outbase", action="store", type="string", dest="outbase", help="output file(s) path/filename base", metavar="FILE")
parser.add_option("-d", "--ssd", action="store", type="string", dest="ssdfile", help="stim separation (SSD) path/filename", metavar="SSD")
parser.add_option("-n", "--cleanoutput", action="store_true", dest="cleanoutput", help="only write unmasked, fully calibrated data -- DEPRECATED")
parser.add_option("-u", "--nullout", action="store_true", dest="nullout", help="write NULL entries to a separate file", metavar="NULL")
parser.add_option("-v", "--verbose", action="store", type="float", dest="verbose", help="Display more output. Set to 0-2", metavar="VRB")
parser.add_option("--retries", action="store", type="int", dest="retries", default=20, help="Set the number of times to ping the server for a response before defining a query failure.  Default is 20, set to a large number if you expect, or want to allow, the query to take a long time.")

(options, args) = parser.parse_args()

#mandatory = ['raw6file', 'aspfile', 'scstfile', 'ssdfile', 'calpath', 'outbase']
mandatory = ['raw6file', 'scstfile', 'calpath', 'outbase']
for m in mandatory:
	if not options.__dict__[m]:
		print "A mandatory option is missing:",m
		parser.print_help()
		os._exit(1)

aspfile = 0
if options.__dict__['aspfile']:
	aspfile = str(options.aspfile).split(',')

ssdfile = 0
if options.__dict__['ssdfile']:
	ssdfile = str(options.ssdfile)

# If the band is not explicity called, attempt to derive it from the raw6 filename.
if not options.band:
	print "Determining band from raw6 filename..."
	band = find_band(options.raw6file)
	if not band:
		print "Unable to parse band from raw6 filename. Specify band on command line using --band."
		parser.print_help()
		os._exit(1)
else:
	band = options.band
	band = band.upper()
	if band != "NUV" and band != "FUV":
		print "Band must be NUV or FUV. "
		os._exit(1)

if not options.cleanoutput:
	options.cleanoutput = 0
if not options.nullout:
	options.nullout = 0
if not options.verbose:
	verbose=0
else:
	verbose=options.verbose

PhotonPipe(options.raw6file,options.scstfile,options.calpath,band,options.outbase,aspfile,ssdfile,options.nullout,verbose=verbose,retries=options.retries)
