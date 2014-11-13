#!/usr/bin/python
import os
import ast
import sys
import argparse
import numpy as np
import curvetools as ct
import gphoton_args as gargs
#from curvetools import
from imagetools import write_jpeg # For JPEG preview image creation
from dbasetools import fGetTimeRanges, suggest_parameters

def gAperture(band,skypos,radius,csvfile=False,annulus=None, coadd=False,
              stepsz=False,verbose=0,clobber=False,trange=None,minexp=1.,
              maxgap=1.,maskdepth=20.,maskradius=1.5,iocode='wb'):
    """Runs gAperture and returns the data in a python dict() and as
    a CSV file if outfile is specified. Can be called from the interpreter.
    """
    if verbose>1:
        print "Generating a light curve with the following paramters:"
        print " band:    {band}".format(band=band)
        print " skypos:  {skypos}".format(skypos=skypos)
        print " trange:  {trange}".format(trange=trange)
        print " radius:  {radius}".format(radius=radius)
        print " annulus: {annulus}".format(annulus=annulus)
        print " stepsz:  {stepsz}".format(stepsz=stepsz)
        print " csvfile: {csvfile}".format(csvfile=csvfile)
        print " verbose: {verbose}".format(verbose=verbose)
    data = ct.write_curve(band, skypos[0], skypos[1], radius, csvfile=csvfile,
                          annulus=annulus, stepsz=stepsz, verbose=verbose,
                          clobber=clobber, trange=trange, coadd=coadd,
                          minexp=minexp, maxgap=maxgap, iocode = iocode,
                          maskdepth=maskdepth, maskradius=maskradius)
    return data

def check_radius(args):
    """Checks the radius value."""
    if not args.radius and not args.suggest:
        print "Must specify an aperture radius."
        exit(0)
    return args

def check_annulus(args):
    """Checks and formats the annulus values."""
    if not (args.annulus1 and args.annulus2) and not args.annulus:
        args.annulus = None
    elif args.annulus1 and args.annulus2:
        args.annulus = [args.annulus1,args.annulus2]
    else:
        args.annulus = list(args.annulus)
    return args

# TODO: fGetTimeRanges is probalby not actually necessary. What we want
# instead is for a way to specifically include or exclude time ranges.
def check_tranges(args):
    """Checks and formats the tranges values."""
    args.tranges = (list(args.tranges) if
                                    args.tranges else [args.tmin,args.tmax])
    args.trange = [np.array(args.tranges).flatten().min(),
                   np.array(args.tranges).flatten().max()]
    if args.verbose:
        print 'Using all exposure in [{tmin},{tmax}]'.format(
                                                tmin=args.tmin,tmax=args.tmax)
    args.tranges = fGetTimeRanges(args.band,args.skypos,maxgap=args.gap,verbose=args.verbose,minexp=args.minexp,trange=args.tranges,detsize=args.detsize)
    if not len(args.tranges):
        print 'No exposure time in database.'
        exit(0)
    if args.verbose:
        print 'Using {t} seconds (raw) in {n} distinct exposures.'.format(
            t=(args.tranges[:,1]-args.tranges[:,0]).sum(),n=len(args.tranges))
    # Make sure tranges is a 2D array
    if len(np.array(args.tranges).shape)==1:
        args.tranges=[args.tranges]
    return args

def setup_parser(iam='gaperture'):
    """Defines command line arguments."""
    parser = argparse.ArgumentParser(description="Generate a light curve")
    parser = gargs.common_args(parser,iam)
    parser.add_argument("-a", "--aperture", action="store", dest="radius",
        help="Aperture radius in decimal degrees",
        type=float)
    parser.add_argument("-i", "--inner", action="store", dest="annulus1",
        help="Inner annulus radius for background subtraction", type=float)
    parser.add_argument("-o", "--outer", action="store", dest="annulus2",
        help="Outer annulus radius for background subtraction", type=float)
    parser.add_argument("--annulus", action="store", type=ast.literal_eval,
        dest="annulus",
        help="Annulus inner and outer radius definition as '[inner,outer]'")
    parser.add_argument("-f", "--file", "--outfile", "--csvfile",
        action="store", type=str, dest="csvfile", help="CSV output file")
    parser.add_argument("--stamp", action="store", type=str, dest="stamp",
        help="Filename for a JPEG preview stamp of the targeted region.")
    parser.add_argument("--addhdr", action="store_true", dest="addhdr",
        help="Add command line and column names to the top of the .csv file.")
    parser.add_argument("--bestparams", "--best", action="store_true",
        dest="best",
        help="Auto set params to produce the highest quality lightcurve.",
        default=False)
    parser.add_argument("--iocode", action="store", dest="iocode", default="wb",
        help="The iocode to be past to the cvs writer. Don't much with this.",
        type=str)
    parser.add_argument("--bgmaskdepth", action="store", dest="maskdepth",
        help="Depth of the background mask in AB Magnitudes.",
        type=float, default=20.0)
    parser.add_argument("--bgmaskradius", action="store", dest="maskradius",
        help="Radius of background mask in n sigmas (assuming Gaussians)",
        type=float, default=1.5)
    return parser

def check_args(args,iam='gaperture'):
    """Checks validity of command line arguments and, in some cases
    massages them a little bit.
    """
    args = gargs.check_common_args(args,iam)
    args = check_radius(args)
    args = check_annulus(args)
    #args = check_tranges(args)
    return args

def reconstruct_command(args):
    """Reconstruct the command line."""
    cmd = " ".join(sys.argv)
    if args.verbose > 1:
        print cmd
    return cmd

def setup_file(args):
    """ If requested, create a header for the CSV that includes the column
    names and a reconstruction of the command line call.
    """
    if not args.csvfile:
        args.csvfile = False
    else:
    # This initiates the file.
        if not args.overwrite and os.path.exists(args.csvfile):
            print "{csvfile} exists. Pass --overwrite to replace it.".format(
                                                          csvfile=args.csvfile)
            exit(0)
        f = open(args.csvfile,args.iocode)
        if args.addhdr:
            if args.verbose:
                print 'Recording command line construction to {f}'.format(
                                                                f=args.csvfile)
            # Write the command line to the top of the output file
            f.write('| '+reconstruct_command(args)+'\n')
        f.close()
        # Setting iocode to append to the file we just created
        args.iocode = 'ab'
    return args

def stamp(args):
    if not args.stamp:
        return
    else:
        if args.verbose:
            print 'Writing stamp preview to {f}'.format(f=args.stamp)
        width = 2.*(args.annulus2 if args.annulus2 else args.radius)
        args.skyrange = [width,width]
        write_jpeg(args.stamp, args.band, args.skypos, args.tranges,
                   args.skyrange, width=False, height=False,
                   stepsz=args.stepsz, clobber=args.overwrite,
                   verbose=args.verbose)
        return

if __name__ == '__main__':
    """Called when gAperture is executed through the command line."""
    args = setup_parser().parse_args()
    args = check_args(args)
    args = setup_file(args)
    stamp(args)
    # TODO: add support for trange(s)
    data = gAperture(args.band, args.skypos, args.radius, csvfile=args.csvfile,
                     annulus=args.annulus, stepsz=args.stepsz,
                     verbose=args.verbose, clobber=args.overwrite,
                     trange=args.trange, coadd=args.coadd, minexp=args.minexp,
                     maxgap=args.maxgap, iocode=args.iocode,
                     maskdepth=args.maskdepth, maskradius=args.maskradius)
