#!/usr/bin/python
import os
import ast
import argparse
from imagetools import * #FIXME: dangerous import
import dbasetools as dbt
import gphoton_args as gargs
import numpy as np

def gMap(band=None,cntfile=False,
		 coadd=False,detsize=1.25,intfile=False,rrfile=False,skypos=None,
		 maxgap=1500.,memlight=100.,minexp=1.,overwrite=False,retries=20,
		 skyrange=None,stepsz=0.,trange=None,verbose=0):
	"""Use a mix of strings (if we want to make an output file) and Booleans
	(False if we do not). I don't think this is the best way (I'd rather
	have separate variables for the True/False to create an image, and a
	string with the name of the output image), but for now this is kept since
	this is how the code was originally written.
	"""
	# FIXME: Maybe?
	response = True if (intfile or rrfile) else False
	write_cnt = cntfile if (cntfile) else False
	write_int = intfile if (intfile) else False
	write_rr = rrfile if (rrfile) else False

	write_images(band,skypos,trange,args.skyrange,width=False,height=False,
				 write_cnt=write_cnt,write_int=write_int,write_rr=write_rr,
				 framesz=stepsz,clobber=overwrite,verbose=verbose,
				 memlight=memlight,coadd=coadd,response=response,
				 retries=retries)

def setup_parser(iam='gmap'):
	parser = argparse.ArgumentParser(description="Generate images / maps.")
	parser = gargs.common_args(parser,iam)
	parser.add_argument("--angle", action="store", type=float, dest="angle",
		default=None, help="The angle subtended in both RA and DEC in degrees.")
	parser.add_argument("--count", action="store", type=str, dest="cntfile",
		default=None, help="File name (full path) for the count image.")
	parser.add_argument("--raangle", action="store", type=float, dest="raangle",
		help="The angle of sky in degrees that the right ascension subtends. "+
		"Overrides --angle.",default=None)
	parser.add_argument("--decangle", action="store", type=float,
		dest="decangle", default=None,
		help="The angle of sky in degrees that the declination subtends. "+
		"Overrides --angle.")
	parser.add_argument("--skyrange", action="store", dest="skyrange",
		type=ast.literal_eval, help="Two element list of ra and dec ranges. "+
		"Equivalent to separately setting --raangle and decangle.")
	parser.add_argument("--intensity", action="store", type=str, dest="intfile",
		default=None, help="File name (full path) for the intensity image.")
	parser.add_argument("--memlight", action="store", type=float,
		dest="memlight", default=100., help="Reduce server-side memory usage "+
		"by requesting data in chunks of no more than this depth in seconds.")
	parser.add_argument("--response", action="store", type=str, dest="rrfile",
		help="File name (full path) for the response image.", default=None)
	return parser

def check_args(args,iam='gmap'):
	args = gargs.check_common_args(args,iam)

	# Check that either skypos _or_ (raangle _and_ decangle) are defined.
	if (not args.skyrange and not (args.raangle and args.decangle)
															and not args.angle):
		raise SystemExit('Must specify either skyrange or ra/dec angle.')
	if (args.raangle and not args.decangle) or (not args.raangle and args.decangle):
		raise SystemExit('Must specify both raangle and deangle or neither.')

	# --angle overwrites everything
	if args.angle:
		args.skyrange = [args.angle,args.angle]

	# --skryange overwrites everything else
	if args.skyrange:
		if np.array(args.skyrange).shape==(2,):
			args.raangle, args.decangle = args.skyrange
		else:
			gPhotonArgsError("Invalid --skyrange: {s}".format(s=args.skyrange))

	# use --angle to fill in missing subtend angles
	if args.angle and not args.raangle:
		args.raangle = args.angle
	if args.angle and not args.decangle:
		args.decangle = args.angle

	# check that requested image subtend angles are positive
	for i,a in enumerate([args.angle,args.raangle,args.decangle]):
		if not a:
			continue
		if a<=0:
			raise SystemExit("Subtended angles must be > 0 degrees.")

	args.skyrange = [args.raangle,args.decangle]

	if args.detsize and args.detsize <= 0.:
		raise SystemExit("Mask radius must be > 0 degrees.")
	if args.memlight and args.memlight <= 0.:
		raise SystemExit("Maximum data chunk per must be > 0 seconds.")

	for image in [args.cntfile, args.intfile, args.rrfile]:
		# Check for overwriting existing images.
		if image and os.path.exists(image) and not args.overwrite:
			raise SystemExit("{f} already exists.".format(f=image))
		# Check if you need to create a new directory and create it
		if image and not os.path.isdir(os.path.dirname(image)):
			print 'Creating directory: {d}'.format(
									d=os.path.abspath(os.path.dirname(image)))
			os.makedirs(os.path.abspath(os.path.dirname(image)),0755)

	return args

if __name__ == '__main__':
	"""Called when gMap is executed directly through command line."""
	args = setup_parser().parse_args()
	args = check_args(args)
	if args.verbose:
		print '{a} image(s) at {pos}'.format(
				a='Writing' if not args.coadd else 'Coadding', pos=args.skypos)
		print '		of dimensions {w}x{h} degrees'.format(w=args.skyrange[0],
														  h=args.skyrange[1])
		print '		in time range(s): {t}'.format(t=args.trange)
	gMap(band=args.band, cntfile=args.cntfile,
		coadd=args.coadd, detsize=args.detsize, intfile=args.intfile,
		rrfile=args.rrfile, skypos=args.skypos, maxgap=args.maxgap,
		memlight=args.memlight, minexp=args.minexp, overwrite=args.overwrite,
		retries=args.retries, skyrange=args.skyrange, stepsz=args.stepsz,
		trange=args.trange, verbose=args.verbose)
