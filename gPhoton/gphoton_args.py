"""
.. module:: gphoton_args
   :synopsis: This module contains functions for checking duplicate arguments
       across gFind, gAperture, gMap, and other gPhoton functions.

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
import ast
from builtins import str
import os
import numpy as np
# gPhoton imports.
import gPhoton.dbasetools as dbt

# ------------------------------------------------------------------------------
class gPhotonArgsError(Exception):
    """
    Exception class specific to gphoton_args.
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def common_args(parser, function_name,
                valid_functions=['gaperture', 'gmap', 'gfind']):
    """
    Defines the arguments and options for the parser object when
        called from the command line.  Accepts a string used to determine which
        arguments to add to the Parser object. Valid function names are "gfind",
        "gaperture", or "gmap" (all case insensitive).

    :param parser: Command-line options to check and modify.

    :type parser: argparse.ArgumentParser Namespace

    :param function_name: Name of the function being called.

    :type function_name: str

    :param valid_functions: List of known/valid functions.

    :type valid_functions: list

    :returns: argparse.ArgumentParser Namespace -- The updated command-line
        arguments.
    """

    try:
        function_name = function_name.strip().lower()
    except AttributeError:
        raise gPhotonArgsError("Invalid function: {f}".format(f=function_name))

    if not function_name in valid_functions:
        raise gPhotonArgsError("{f} not in {vf}".format(f=function_name,
                                                        vf=valid_functions))

    parser.add_argument(
        "-b", "--band", action="store", type=lambda b: str.upper(str(b)),
        dest="band", help="Band designation",
        default=str(u"BOTH") if function_name == 'gfind' else str(u"NUV"),
        choices=[str(u"NUV"), str(u"FUV")]+(
            [str(u"BOTH")] if function_name == 'gfind' else []))

    parser.add_argument("-d", "--dec", action="store", type=float,
                        dest="dec", metavar="DEC",
                        help="Center Declination position in decimal degrees. "
                        "Must be 0 < DEC < 90.")
    parser.add_argument("--detsize", action="store", type=float,
                        dest="detsize", default=1.1,
                        help="Set the effective field diameter in degrees for"
                        " the exposure search.  Default = 1.1.")
    parser.add_argument("-g", "--gap", "--maxgap", action="store",
                        type=float, dest="maxgap", default=1500.,
                        help="Maximum gap size in seconds for data to be"
                        " considered contiguous.  Default = 1500.")
    parser.add_argument("--minexp", action="store", type=float,
                        dest="minexp", help="Minimum contiguous exposure in"
                        " seconds for data to be reported.  Default = 1.",
                        default=1.)
    parser.add_argument("-r", "--ra", action="store", type=float, dest="ra",
                        help="Center Right Ascension position in decimal"
                        " degrees. Must be 0 < RA < 360.", metavar="RA")
    parser.add_argument("--retries", action="store", type=int, dest="retries",
                        help="Set the number of times to ping the server for a"
                        " response before defining a query failure. Default is"
                        " 20, set to a large number if you expect, or want to"
                        " allow, the query to take a long time without issuing"
                        " a failure.", default=20)
    parser.add_argument("--skypos", action="store", dest="skypos",
                        help="Alternate method for specifying sky position"
                        " with format '[RA,Dec]'", type=ast.literal_eval)
    parser.add_argument("--t0", "--tmin", action="store", type=float,
                        dest="tmin", help="Minimum date of observation to"
                        " consider (specify in GALEX time standard).  Default"
                        " = 6e8", default=6.e8)
    parser.add_argument("--t1", "--tmax", action="store", type=float,
                        dest="tmax", help="Maxium date of observation to"
                        " consider (specify in GALEX time standard).  Default"
                        " = 11e8", default=11.e8)
    parser.add_argument("--trange", "--tranges", action="store", dest="trange",
                        help="Time range(s) in which to limit the search, in"
                        " the format '[t0,t1]' or '[[t0_a,t1_a],[t0_b,t1_b]]'"
                        " (format in GALEX time).", type=ast.literal_eval)
    parser.add_argument("-v", "--verbose", action="store", type=int,
                        dest="verbose", help="Prints extra information to"
                        " STDOUT (higher number = more output). Choices are"
                        " {0,1,2,3}, default = 0.", default=0,
                        choices=[0, 1, 2, 3])
    parser.add_argument("--suggest", action="store_true", dest="suggest",
                        help="Suggest reasonable parameters for aperture"
                        " photometry. The includes recenting on the nearest"
                        " MCAT source. This flag will clobber other annuli and"
                        " aperture radii parameters.", default=False)
    parser.add_argument("--skyrange", action="store", dest="skyrange",
                        type=ast.literal_eval, help="Two element list of ra"
                        " and dec ranges. Equivalent to separately setting"
                        " --raangle and decangle.")

    if function_name in ['gaperture', 'gmap']:
        parser.add_argument("--calpath", action="store", type=str,
                            dest="calpath",
                            default=os.pardir+os.sep+"cal"+os.sep,
                            help="Path to the directory that contains the"
                            " calibration files.")
        parser.add_argument("--coadd", action="store_true", dest="coadd",
                            help="Return the coadded flux (gAperture) or a"
                            " coadded image (gMap) over all requested time"
                            " ranges?  Default = False.", default=False)
        parser.add_argument("--overwrite", "--ow", "--clobber",
                            action="store_true", dest="overwrite",
                            help="Overwrite existing output files?  Default ="
                            " False.", default=False)
        parser.add_argument("-s", "--step", "--stepsz", "--frame",
                            action="store", type=float, dest="stepsz",
                            help="Step size for lightcurve or movie in"
                            " seconds.  Default = 0. (no binning).", default=0.)

    return parser
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def check_common_args(args, function_name,
                      valid_functions=['gaperture', 'gmap', 'gfind'],
                      allow_no_coords=False):
    """
    Checks validity of some command line arguments used in gFind,
        gAperture, gMap, etc.  Returns the appropriate arguments as variables
        back to the calling procedure.

    :param args: Command-line options to check and modify.

    :type args: argparse.ArgumentParser Namespace

    :param function_name: Name of the function being called.

    :type function_name: str

    :param valid_functions: List of known/valid functions.

    :type valid_functions: list

    :param allow_no_coords: Allow function to be called without coordinates?

    :type allow_no_coords: bool

    :returns: argparse.ArgumentParser Namespace -- The updated command-line
        arguments.
    """

    try:
        function_name = function_name.strip().lower()
    except AttributeError:
        raise gPhotonArgsError("Invalid function: {f}".format(f=function_name))

    if not function_name in valid_functions:
        raise gPhotonArgsError("Invalid function: {f}".format(f=function_name))

    try:
        args.band = args.band.strip()
    except AttributeError:
        raise SystemExit("Invalid band: {b}".format(b=args.band))

    # This will ensure calpath has a trailing '/'.
    if function_name in ['gaperture', 'gmap']:
        args.calpath = os.path.join(args.calpath, '')
        # [Future]: Consider fixing this statement. This is breaking nosetests,
        # but it's not a bad idea...
        # if not os.path.isdir(args.calpath):
        #    raise SystemExit("Calibration path not found: " + args.calpath)

    if (not (args.ra and args.dec) and not args.skypos and
            not allow_no_coords):
        raise SystemExit("Must specify either both RA/DEC or SKYPOS.")
    elif (args.ra and args.dec) and args.skypos:
        if not (args.ra == args.skypos[0] and args.dec == args.skypos[1]):
            raise SystemExit("Must specify either RA/DEC or SKYPOS, not both.")
    elif (args.ra and args.dec) and not args.skypos:
        args.skypos = [args.ra, args.dec]
    elif not (args.ra and args.dec) and args.skypos:
        args.ra, args.dec = args.skypos

    if args.suggest and function_name in ['gfind', 'gaperture']:
        (args.ra, args.dec, args.radius, args.annulus1,
         args.annulus2) = dbt.suggest_parameters(args.band, args.skypos,
                                                 verbose=0)
        args.skypos = [args.ra, args.dec]
        if args.verbose:
            print("Recentering on ["+str(args.ra)+", "+str(args.dec)+"]")
            print("Setting radius to "+str(args.radius))
            print("Setting annulus to ["+str(args.annulus1)+", "+
                  str(args.annulus2)+"]")

    if args.skypos:
        if np.array(args.skypos).shape != (2,):
            raise gPhotonArgsError(
                "Skypos (--skypos) must be a 2-element array.")
        args.ra, args.dec = args.skypos

    if args.ra and not 0. <= args.ra <= 360.:
        raise SystemExit(
            "RA of {ra} does not satisfy 0 <= RA <= 360".format(ra=args.ra))

    if args.dec and not -90 <= args.dec <= 90:
        raise SystemExit(
            "Dec of {dec} does not satisfy -90 <= DEC <= 90".format(
                dec=args.dec))

    if args.detsize and args.detsize <= 0.:
        raise SystemExit("Effective field diameter (--detsize) must be > 0")

    if args.maxgap and args.maxgap <= 0.:
        raise SystemExit("Maximum gap length (--maxgap) must be > 0 seconds.")
    if args.minexp and args.minexp <= 0.:
        raise SystemExit("Minimum valid exposure depth (--minexp) must be > 0"
                         " seconds.")

    if args.retries and args.retries <= 0.:
        raise SystemExit("Number of retries (--retries) must be > 0.")

    # tmin / tmax must be defined and reasonable
    if not args.tmin or args.tmin <= 0.:
        raise SystemExit("T0 (--t0) must be > 0.")
    if not args.tmax or args.tmax <= 0.:
        raise SystemExit("T1 (--t1) must be > 0.")
    if args.tmin >= args.tmax:
        raise SystemExit("Minimum time (--t0) must be < maximum time (--t1).")

    if args.trange:
        if np.array(args.trange).shape == (2, ):
            args.trange = [args.trange]
        if not (len(np.array(args.trange).shape) == 2 and
                np.array(args.trange).shape[1] == 2):
            raise SystemExit("trange (--trange) must be a pairwise list.")
        # Individually check the entries for sanity
        for t in args.trange:
            if t[0] <= 0 or t[1] <= 0:
                raise SystemExit('Times must be positive: {t}'.format(t=t))
            if t[1] <= t[0]:
                raise SystemExit('Start time ({t0}) must preceed end time'
                                 ' ({t1})'.format(t0=t[0], t1=t[1]))
    elif not allow_no_coords and function_name in ['gmap', 'gaperture']:
        args.trange = dbt.fGetTimeRanges(args.band, args.skypos,
                                         trange=[args.tmin, args.tmax],
                                         maxgap=args.maxgap, minexp=args.minexp,
                                         detsize=args.detsize,
                                         skyrange=args.skyrange)
    else:
        # If no coordinates specified then use a huge time range for now.
        args.trange = [args.tmin, args.tmax]

    return args
# ------------------------------------------------------------------------------
