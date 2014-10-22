#!/usr/bin/python
import os
import ast
import argparse
import numpy as np
import curvetools as ct
#from curvetools import
from imagetools import write_jpeg # For JPEG preview image creation
from optparse import OptionParser
from dbasetools import fGetTimeRanges, suggest_parameters
import sys

#def gAperture(args):
def gAperture(band,skypos,radius,csvfile=False,annulus=None, coadd=False,
              stepsz=False,verbose=0,clobber=False,trange=None,minexp=1.,
              maxgap=1.):
    """Runs gAperture and returns the data in a python dict() and as
    a CSV file if outfile is specified. Can be called from the interpreter.
    """
    if verbose>1:
        print "Generating a light curve with the following paramters:"
        print " band:    {bana}".format(band=band)
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
                          minexp=minexp, maxgap=maxgap)
    return data

def check_radius(args):
    """Checks the radius value."""
    if not args.radius and not args.suggest:
        print "Must specify an aperture radius."
        exit(0)
    return

def suggest(args):
    """Generates suggested lightcurve parameters with suggest_parameters()"""
    (args.ra, args.dec, args.radius, args.annulus1,
                   args.annulus2) = suggest_parameters(args.band,args.skypos)
    if args.verbose:
        print "Recentering on [{ra},{dec}]".format(ra=args.ra,dec=args.dec)
        print "Setting radius to {radius}".format(radius=args.radius)
        print "Setting annulus to [{ann1},{ann2}]".format(
                                        ann1=args.annulus1,ann2=args.annulus2)
    return args

def check_skypos(args):
    """Checks and formats the skypos values."""
    if not (args.ra and args.dec) and not args.skypos:
        print "Must specify either RA/Dec or skypos."
        exit(0)
    elif (args.ra and args.dec) and args.skypos:
        print "Must specify either RA/Dec or skypos. Not both."
        exit(0)
    elif (args.ra and args.dec):
        args.skypos = [args.ra,args.dec]
    else:
        args.skypos = list(args.skypos)
    if args.suggest:
        args = suggest(args)
    return args

def check_stepsz(args):
    """Checks the stepsz value."""
    if args.stepsz and args.coadd:
        print "Cannot specify both --stepsz and --coadd."
        exit(0)
    return

def check_annulus(args):
    """Checks and formats the annulus values."""
    if not (args.annulus1 and args.annulus2) and not args.annulus:
        args.annulus = False
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

def check_args(args):
    """Checks validity of command line arguments and, in some cases
    massages them a little bit.
    """
    check_radius(args)
    check_stepsz(args)
    args = check_skypos(args)
    args = check_annulus(args)
    args = check_tranges(args)
    return args

def setup_parser():
    """Defines command line arguments."""
    parser = argparse.ArgumentParser(description="Generate a light curve")
    parser.add_argument("-b", "--band", action="store", dest="band",
        help="[NF]UV band designation", choices=['NUV','FUV'], type=str.upper,
        required=True)
    parser.add_argument("-r", "--ra", action="store", dest="ra",
        help="Center Right Ascension position in decimal degrees",
        type=np.float64)
    parser.add_argument("-d", "--dec", action="store", dest="dec",
        help="Center Declination position in decimal degrees", type=np.float64)
    parser.add_argument("--skypos", action="store", dest="skypos",
        help="Alternate method for specifying sky position as '[RA,Dec]'",
        type=ast.literal_eval)
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
    parser.add_argument("--t0", "--tmin", action="store", dest="tmin",
        help="Minimum time to consider.",default=1.,type=np.float64)
    parser.add_argument("--t1", "--tmax", action="store", dest="tmax",
        help="Maxium time to consider.",default=1000000000000.,type=np.float64)
    parser.add_argument("--tranges", "--trange", action="store",
        type=ast.literal_eval, dest="tranges",
        help="List of time ranges with format '[[t0,t1],[t2,t3],...]'")
    parser.add_argument("-s", "--step", "--frame", action="store", type=float,
        dest="stepsz", help="Step size in seconds",default=0)
    parser.add_argument("-f", "--file", "--outfile", "--csvfile",
        action="store", type=str, dest="csvfile", help="CSV output file")
    parser.add_argument("-v", "--verbose", action="store", type=int,
        dest="verbose", help="Display more output. Set to 0-2.", default=0,
        choices=[0,1,2,3])
    parser.add_argument("--stamp", action="store", type=str, dest="stamp",
        help="Filename for a JPEG preview stamp of the targeted region.")
    parser.add_argument("--calpath", action="store", type=str, dest="calpath",
        help="Path to the directory that contains the calibration files.",
        default='../cal/')
    parser.add_argument("--addhdr", action="store_true", dest="addhdr",
        help="Add command line and column names to the top of the .csv file.")
    parser.add_argument("--coadd", action="store_true", dest="coadd",
        help="Return the coadded flux over all requested time ranges.")
    parser.add_argument("-g", "--gap", "--maxgap", action="store",
        type=float, dest="gap",
        help="Max gap size in seconds to be considered contiguous.",
        default=1)
    parser.add_argument("--minexp", action="store", type=float,
        dest="minexp",
        help="Minimum contiguous exposure in seconds for data to be reported.",
        default=1)
    parser.add_argument("--detsize", action="store", type=float,
        dest="detsize",
        help="Set the FOVdiameter in degrees for the exposure search.",
        default=1.25)
    parser.add_argument("--bestparams", "--best", action="store_true",
        dest="best",
        help="Auto set params to produce the highest quality lightcurve.",
        default=False)
    parser.add_argument("--suggest", "--optimize", action="store_true",
        dest="suggest",
        help="Suggest optimum parameters for aperture photometry.",
        default=False)
    parser.add_argument("--overwrite", "--ow", "--clobber",
        action="store_true", dest="overwrite", default=False,
        help="Overwrite any preexisting files. Will supress warnings.")
    parser.add_argument("--iocode", action="store", dest="iocode",
        default="wb",
        help="The iocode to be past to the cvs writer. Don't much with this.",
        type=str)
    return parser

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
            print "File {csvfile} exists. Pass --overwrite to replace it.".format(csvfile=args.csvfile)
            exit(0)
        f = open(args.csvfile,args.iocode)
        if args.addhdr:
            # Write the command line to the outfile; should be optional
            f.write('| '+reconstruct_command(args)+'\n')
            # FIXME: These are no longer the correct column names.
            f.write('| tstart, tstop, radius, exptime, cps, error, flux, flux_error, magnitude, mag_error, inner annulus, outer annulus, background, response, counts, aperture correction 1, aperture correction 2\n')
        f.close()
        # Setting this will now append to the file we just created
        args.iocode = 'ab'
    return args

def stamp(args):
    if not args.stamp:
        return
    else:
        if args.verbose > 1:
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
                     maxgap=args.gap)
