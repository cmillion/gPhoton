#!/usr/bin/python

import os

from PhotonPipe import *

from optparse import OptionParser

def gPhoton(raw6file,scstfile,band,outbase,aspfile,ssdfile,nullout,
			verbose=0,retries=20):
	PhotonPipe(raw6file,scstfile,band,outbase,aspfile,ssdfile,nullout,
			   verbose=verbose,retries=retries)
	return

def setup_parser():
	parser = OptionParser()
	parser.add_option("-r", "--raw6", action="store", type="string",
		dest="raw6file", help="raw6 path/filename", metavar="FILE")
	parser.add_option("-a", "--aspect", action="store", type="string",
		dest="aspfile", help="aspection path/filename", metavar="FILE")
	parser.add_option("-s", "--scst", action="store", type="string",
		dest="scstfile", help="spacecraft state path/filename", metavar="SCST")
	parser.add_option("-b", "--band", action="store", type="string",
		dest="band", help="[NF]UV band designation", metavar="BAND")
	parser.add_option("-o", "--outbase", action="store", type="string",
		dest="outbase", help="output file(s) path/filename base",
		metavar="FILE")
	parser.add_option("-d", "--ssd", action="store", type="string",
		dest="ssdfile", help="stim separation (SSD) path/filename",
		metavar="SSD")
	parser.add_option("-n", "--cleanoutput", action="store_true",
		dest="cleanoutput", default=0,
		help="only write unmasked, fully calibrated data -- DEPRECATED")
	parser.add_option("-u", "--nullout", action="store_true", dest="nullout",
		help="write NULL entries to a separate file", metavar="NULL",
		default=0)
	parser.add_option("-v", "--verbose", action="store", type="float",
		dest="verbose", help="Display more output. Set to 0-2", metavar="VRB",
		default=0)
	parser.add_option("--retries", action="store", type="int", dest="retries",
		default=20, help="Query attempts before timeout.")
	return parser

def check_args(args):
	try:
		if not (args.raw6file and args.scstfile and args.outbase):
			print "A mandatory option is missing."
			os._exit(1)
	except:
		print "A mandatory option is missing: raw6file, scstfile, or outbase"
		os._exit(1)

	if args.aspfile:
		args.aspfile = str(args.aspfile).split(',')

	if args.ssdfile:
		args.ssdfile = str(args.ssdfile)

	# If the band is not explicity called, attempt to derive it from the raw6 filename.
	if not args.band:
		print "Determining band from raw6 filename..."
		args.band = find_band(args.raw6file)
		if not args.band:
			print "Unable to parse band from raw6 filename. Specify band on command line using --band."
			os._exit(1)
	else:
		args.band = args.band.upper()
		if not band in ["NUV","FUV"]:
			print "Band must be NUV or FUV. "
			os._exit(1)

	return args

def __main__():
	args = setup_parser().parse_args()
	args = check_args(args)
	gPhoton(args.raw6file,args.scstfile,args.band,args.outbase,args.aspfile,
		args.ssdfile,args.nullout,verbose=args.verbose,retries=args.retries)

if __name__ == "__main__":
    try:
        __main__()
    except (KeyboardInterrupt, pycurl.error):
        exit('Received Ctrl + C... Exiting! Bye.', 1)
