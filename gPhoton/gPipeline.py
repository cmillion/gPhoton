#!/usr/bin/python

import os
from PhotonPipe import PhotonPipe
import argparse

def gPipeline(raw6file,scstfile,band,outbase,aspfile,ssdfile,nullout,
			verbose=0,retries=20):
	PhotonPipe(raw6file,scstfile,band,outbase,aspfile,ssdfile,nullout,
			   verbose=verbose,retries=retries)
	return

def setup_parser():
	parser = argparse.ArgumentParser(description="Generate time-tagged lists of aspect-corrected photons positions.")
	parser.add_argument("-r", "--raw6", action="store", type=str,
		dest="raw6file", help="raw6 path/filename")
	parser.add_argument("-a", "--aspect", action="store", type=str,
		dest="aspfile", help="aspection path/filename")
	parser.add_argument("-s", "--scst", action="store", type=str,
		dest="scstfile", help="spacecraft state path/filename")
	parser.add_argument("-b", "--band", action="store", type=str,
		dest="band", help="[NF]UV band designation")
	parser.add_argument("-o", "--outbase", action="store", type=str,
		dest="outbase", help="output file(s) path/filename base", default=None)
	parser.add_argument("-d", "--ssd", action="store", type=str,
		dest="ssdfile", help="stim separation (SSD) path/filename")
	parser.add_argument("-u", "--nullout", action="store_true", dest="nullout",
		help="write NULL entries to a separate file", default=False)
	parser.add_argument("-v", "--verbose", action="store", type=float,
		dest="verbose", help="Display more output. Set to 0-2", default=0)
	parser.add_argument("--retries", action="store", type=int, default=20,
		help="Query attempts before timeout.")
	return parser

def check_args(args):
	if not args.raw6file:
		print "Must specify a RAW6 filename (--raw6)."
		raise SystemExit
	if not args.scstfile:
		print "Must specify a SCST filename (--scst)."
		raise SystemExit
	if not args.outbase:
		print "Must specify an output base filename (--outbase)."
		raise SystemExit

	if args.aspfile:
		args.aspfile = str(args.aspfile).split(',')

	if args.ssdfile:
		args.ssdfile = str(args.ssdfile)

	# If the band is not explicity called, attempt to derive it from the raw6 filename.
	if not args.band:
		print "Determining band from raw6 filename..."
		if '-fd-raw6' in args.raw6file:
			args.band = 'FUV'
		elif '-nd-raw6' in args.raw6file:
			args.band = 'NUV'
		else:
			print "Unable to parse band from raw6 filename. Specify band on command line using --band."
			raise SystemExit
	else:
		args.band = args.band.upper()
		if not band in ["NUV","FUV"]:
			print "Band must be NUV or FUV. "
			raise SystemExit

	return args

def __main__():
	args = setup_parser().parse_args()
	print args
	args = check_args(args)
	gPipeline(args.raw6file,args.scstfile,args.band,args.outbase,args.aspfile,
		args.ssdfile,args.nullout,verbose=args.verbose,retries=args.retries)

if __name__ == "__main__":
    try:
        __main__()
    except (KeyboardInterrupt, pycurl.error):
        exit('Received Ctrl + C... Exiting! Bye.', 1)
