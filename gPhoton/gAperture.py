#!/usr/bin/python

"""
.. module:: gAperture
   :synopsis: Module for the creation of GALEX  light curves with user-defined
       time bins and photometric apertures.
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
import argparse
import ast
import os
import sys
# gPhoton imports.
import gPhoton.curvetools as ct
from gPhoton.galextools import aper2deg
from gPhoton import __version__
import gPhoton.gphoton_args as gargs
from gPhoton.imagetools import write_jpeg # For JPEG preview image creation

# ------------------------------------------------------------------------------
def gaperture(band, skypos, radius, csvfile=None, annulus=None, coadd=False,
              stepsz=False, verbose=0, overwrite=False, trange=None,
              tranges=None, minexp=1., maxgap=1500., iocode='w',
              detsize=1.1, minimal_output=False, photoncsvfile=None,
              addhdr=False):
    """
    Creates a light curve and returns the data in a python dict() and as
        a CSV file, if outfile is specified. Can be called from the interpreter.

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param radius: The radius of the photometric aperture, in degrees.

    :type radius: float

    :param csvfile: Name of the photon event CSV file to use for lightcurve.

    :type csvfile: str

    :param annulus: Radii of the inner and outer of an annulus, in degrees,
        within which to measure the background.

    :type annulus: list

    :param coadd: Set to True if calculating a total flux instead of flux
        from each time bin.

    :type coadd: bool

    :param stepsz: The size of the time bins to use, in seconds.

    :type stepsz: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param overwrite: If True, overwite an existing output file.

    :type overwrite: bool

    :param trange: Minimum and maximum time range to make a light curve,
        in GALEX time seconds.

    :type trange: list

    :param tranges: Set of time ranges to query within, in GALEX time seconds.

    :type tranges: list

    :param minexp: Minimum gap size, in seconds, for data to be considered
        contiguous.

    :type minexp: float

    :param maxgap: Maximum gap size, in seconds, for data to be considered
        contiguous.

    :type maxgap: float

    :param iocode: The code to use when writing the output file.

    :type iocode: str

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :param minimal_output: If True, produce an output file with a minimum
        number of columns.

    :type minimal_output: bool

    :param photoncsvfile: Name of the photon event CSV file to use for photons.

    :type photoncsvfile: str

    :returns: dict -- The light curve, including input parameters.
    """

    if verbose > 1:
        print("Using v{v} of gAperture.".format(v=__version__))
        print("Generating a light curve with the following paramters:")
        print(" band:    {band}".format(band=band))
        print(" skypos:  {skypos}".format(skypos=skypos))
        print(" tranges: {trange}".format(trange=trange if trange else tranges))
        print(" radius:  {radius}".format(radius=radius))
        print(" annulus: {annulus}".format(annulus=annulus))
        print(" stepsz:  {stepsz}".format(stepsz=stepsz))
        print(" csvfile: {csvfile}".format(csvfile=csvfile))
        print(" verbose: {verbose}".format(verbose=verbose))

    data = ct.write_curve(band.upper(), skypos[0], skypos[1], radius,
                          csvfile=csvfile, annulus=annulus, stepsz=stepsz,
                          verbose=verbose, overwrite=overwrite, trange=trange,
                          tranges=tranges, coadd=coadd, minexp=minexp,
                          maxgap=maxgap, iocode=iocode, detsize=detsize,
                          minimal_output=minimal_output,
                          photoncsvfile=photoncsvfile,addhdr=addhdr)

    return data
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def check_radius(args):
    """
    Checks the radius value.

    :param args: The command-line arguments.

    :type args: argparse.ArgumentParser Namespace

    :returns: argparse.ArgumentParser Namespace -- The updated command-line
        arguments.
    """

    if not (args.radius or args.suggest or args.aperradius):
        print("Must specify an aperture radius.")
        raise SystemExit

    if args.radius and args.aperradius:
        print("Must not specify both --aperture and --mcataper.")
        raise SystemExit

    if args.aperradius and not args.radius:
        args.radius = aper2deg(args.aperradius)

    return args
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def check_annulus(args):
    """
    Checks and formats the annulus values.

    :param args: The command-line arguments.

    :type args: argparse.ArgumentParser Namespace

    :returns: argparse.ArgumentParser Namespace -- The updated command-line
        arguments.
    """

    if not (args.annulus1 and args.annulus2) and not args.annulus:
        args.annulus = None
    elif args.annulus1 and args.annulus2:
        args.annulus = [args.annulus1, args.annulus2]
    else:
        args.annulus = list(args.annulus)
        if not (args.annulus1 and args.annulus2):
            args.annulus1 = args.annulus[0]
            args.annulus2 = args.annulus[1]

    return args
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def setup_parser(iam='gaperture'):
    """
    Defines command line arguments.

    :param iam: The name of the method being called.

    :type iam: str

    :returns: argparse.ArgumentParser Namespace -- The command-line arguments.
    """

    parser = argparse.ArgumentParser(description="Generate a light curve")
    parser = gargs.common_args(parser, iam)
    parser.add_argument("-a", "--aperture", action="store", dest="radius",
                        help="Aperture radius in decimal degrees", type=float)
    parser.add_argument("--mcataper", "--apercode", action="store",
                        dest="aperradius", type=int)
    parser.add_argument("-i", "--inner", action="store", dest="annulus1",
                        help="Inner annulus radius for background subtraction",
                        type=float)
    parser.add_argument("-o", "--outer", action="store", dest="annulus2",
                        help="Outer annulus radius for background subtraction",
                        type=float)
    parser.add_argument("--annulus", action="store", type=ast.literal_eval,
                        dest="annulus", help="Annulus inner and outer radius"
                        " definition as '[inner,outer]'")
    parser.add_argument("-f", "--file", "--outfile", "--csvfile",
                        action="store", type=str, dest="csvfile",
                        help="CSV output file")
    parser.add_argument("--photoncsvfile",
                        action="store", type=str, dest="photoncsvfile",
                        help="CSV output file for photon list")
    parser.add_argument("--stamp", action="store", type=str, dest="stamp",
                        help="Filename for a JPEG preview stamp of the"
                        " targeted region.")
    parser.add_argument("--addhdr", action="store_true", dest="addhdr",
                        help="Add command line and column names to the top of"
                        " the .csv file.")
    parser.add_argument("--iocode", action="store", dest="iocode", default="w",
                        help="The iocode to be passed to the cvs writer."
                        " Don't mess with this.", type=str)
    parser.add_argument("--bgmaskdepth", action="store", dest="maskdepth",
                        help="DEPRECATED. Depth of the background mask in AB"
                        " Magnitudes.", type=float, default=20.0)
    parser.add_argument("--bgmaskradius", action="store", dest="maskradius",
                        help="DEPRECATED. Radius of background mask in n"
                        " sigmas (assuming Gaussians)", type=float, default=1.5)
    parser.add_argument("--minimal_output", "--minout", action="store_true",
                        dest="minimal_output",
                        help=("If csvfile is also set, writes only a small"
                              " number of human-readable columns to the"
                              " lightcurve file."))

    return parser
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def check_args(args, iam='gaperture'):
    """
    Checks validity of command line arguments and, in some cases
        modifies them a little bit.

    :param args: The command-line arguments.

    :type args: argparse.ArgumentParser Namespace

    :param iam: The name of the method being called.

    :type iam: str

    :returns: argparse.ArgumentParser Namespace -- The updated command-line
        arguments.
    """

    args = gargs.check_common_args(args, iam)
    args = check_radius(args)
    args = check_annulus(args)

    return args
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def reconstruct_command(args):
    """
    Rebuild the equivalent command line for this call to gAperture.

    :param args: The command-line arguments.

    :type args: argparse.ArgumentParser Namespace

    :returns: str -- The equivalent command-line for this call to gAperture.
    """

    cmd = " ".join(sys.argv)

    if args.suggest:
        cmd = "{cmd} --annulus [{i},{o}] --aperture {a}".format(
            cmd=cmd, i=args.annulus1, o=args.annulus2, a=args.radius)

    if args.verbose > 1:
        print(cmd)

    return cmd
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def setup_file(args):
    """
    If requested, create a header for the CSV that includes the column
        names and a reconstruction of the command line call.

    :param args: The command-line arguments.

    :type args: argparse.ArgumentParser Namespace

    :returns: argparse.ArgumentParser Namespace -- The updated command-line
        arguments.
    """

    if not args.csvfile:
        args.csvfile = False
    else:
        # This initiates the file.
        if not args.overwrite and os.path.exists(args.csvfile):
            print("{csvfile} exists. Pass --overwrite to replace it.".format(
                csvfile=args.csvfile))
            raise SystemExit
        with open(args.csvfile, args.iocode) as f:
            if args.addhdr:
                if args.verbose:
                    print('Recording command line construction to {f}'.format(
                        f=args.csvfile))
                # Write the command line to the top of the output file
                f.write('| '+reconstruct_command(args)+'\n')
        # Setting iocode to append to the file we just created
        args.iocode = 'a'

    return args
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def stamp(args):
    """
    Creates a jpeg preview image stamp of the targeted region.

    :param args: The command-line arguments.

    :type args: argparse.ArgumentParser Namespace
    """

    if not args.stamp:
        return
    else:
        if args.verbose:
            print('Writing stamp preview to {f}'.format(f=args.stamp))

        width = 2.*(args.annulus2 if args.annulus2 else args.radius)

        args.skyrange = [width, width]

        write_jpeg(args.stamp, args.band, args.skypos, args.trange,
                   args.skyrange, width=False, height=False,
                   stepsz=args.stepsz, overwrite=args.overwrite,
                   verbose=args.verbose)

        return
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def __main__():
    """
    Called when gAperture is executed through the command line.
    """

    args = setup_parser().parse_args()
    args = check_args(args)
    args = setup_file(args)
    stamp(args)

    # [Future]: add support for trange(s).
    data = gaperture(args.band, args.skypos, args.radius, csvfile=args.csvfile,
                     annulus=args.annulus, stepsz=args.stepsz,
                     verbose=args.verbose, overwrite=args.overwrite,
                     trange=[args.tmin, args.tmax], tranges=args.trange,
                     coadd=args.coadd, minexp=args.minexp, maxgap=args.maxgap,
                     iocode=args.iocode, detsize=args.detsize,
                     minimal_output=args.minimal_output,
                     photoncsvfile=args.photoncsvfile,addhdr=args.addhdr)
# ------------------------------------------------------------------------------

if __name__ == "__main__":
    try:
        __main__()
    except KeyboardInterrupt:
        exit('Received Ctrl + C... Exiting! Bye.', 1)
# ------------------------------------------------------------------------------
