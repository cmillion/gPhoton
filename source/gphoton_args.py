""" This module contains functions for checking duplicate arguments across
gFind, gAperture, gMap, and other gPhoton functions.
"""

import argparse
import os
import ast
import dbasetools as dbt
import numpy as np

class gPhotonArgsError(Exception):
    """Exception class specific to gphoton_args."""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def common_args(parser,function_name,
                valid_functions=['gaperture','gmap','gfind']):
    """Defines the arguments and options for the parser object when
    called from the command line.  Accepts a string used to determine which
    arguments to add to the Parser object. Valid function names are "gfind",
    "gaperture", or "gmap" (all case insensitive).
    """
    try:
        function_name = function_name.strip().lower()
    except AttributeError:
        raise gPhotonArgsError("Invalid function name specified")

    if not function_name in valid_functions:
        raise gPhotonArgsError("{f} not in {vf}".format(
                                        f=function_name,vf=valid_functions))

    parser.add_argument("-b", "--band", action="store", type=str,
        dest="band", help="Band designation", default="BOTH" if
        function_name=='gfind' else "NUV",
        choices=["NUV","FUV", "nuv", "fuv", "both"]+(['BOTH','both']
        if function_name=='gfind' else []))

    parser.add_argument("-d", "--dec", action="store", type=float,
        dest="dec", metavar="DEC",
        help="Center Declination position in decimal degrees. "
             "Must be 0 < DEC < 90.")
    parser.add_argument("--detsize", action="store", type=float,
        dest="detsize", default=1.25,
        help="Set the effective field diameter in degrees for the "
             "exposure search.  Default = 1.25.")
    parser.add_argument("-g", "--gap", "--maxgap", action="store",
        type=float, dest="maxgap", default=1500.,
        help="Maximum gap size in seconds for data to be considered "
             "contiguous.  Default = 1500.")
    parser.add_argument("--minexp", action="store", type=float,
        dest="minexp", help="Minimum contiguous exposure in seconds for "+
        "data to be reported.  Default = 1.", default=1.)
    parser.add_argument("-r", "--ra", action="store", type=float, dest="ra",
        help="Center Right Ascension position in decimal degrees. "+
        "Must be 0 < RA < 360.", metavar="RA")
    parser.add_argument("--retries", action="store", type=int, dest="retries",
        help="Set the number of times to ping the server for a response "+
        "before defining a query failure. Default is 20, set to a large "+
        "number if you expect, or want to allow, the query to take a long "+
        "time without issuing a failure.", default=20)
    parser.add_argument("--skypos", action="store", dest="skypos",
        help="Alternate method for specifying sky position with format "+
        "'[RA,Dec]'",type=ast.literal_eval)
    parser.add_argument("--t0", "--tmin", action="store", type=float,
        dest="tmin", help="Minimum date of observation to consider (specify "+
        "in GALEX time standard).  Default = 1.",default=1.)
    parser.add_argument("--t1", "--tmax", action="store", type=float,
        dest="tmax", help="Maxium date of observation to consider (specify "+
        "in GALEX time standard).  Default = 1000000000000.",
        default=1000000000000.)
    parser.add_argument("--trange", "--tranges", action="store", dest="trange",
        help="Time range in which to limit the search, in the format "+
        "'[t0,t1]' (specify in GALEX time standard).", type=ast.literal_eval)
    parser.add_argument("-v", "--verbose", action="store", type=int,
        dest="verbose", help="Prints extra information to STDOUT (higher "+
        "number = more output). Choices are {0,1,2,3}, default = 0.",
        default=0, choices=[0,1,2,3])

    if function_name in ['gaperture','gmap']:
        parser.add_argument("--calpath", action="store", type=str,
            dest="calpath", default=os.pardir+os.sep+"cal"+os.sep,
            help="Path to the directory that contains the calibration files.")
        parser.add_argument("--coadd", action="store_true", dest="coadd",
            help="Return the coadded flux (gAperture) or a coadded image (gMap)"
                 " over all requested time ranges?  Default = False.",
                 default=False)
        parser.add_argument("--overwrite", "--ow", "--clobber",
            action="store_true", dest="overwrite", help="Overwrite existing "+
            "output files?  Default = False.", default=False)
        parser.add_argument("-s", "--step", "--frame", action="store",
            type=float, dest="stepsz", help="Step size for lightcurve or "+
            "movie in seconds.  Default = 0. (no binning).", default=0.)

    if function_name in ['gfind','gaperture']:
        parser.add_argument("--suggest", action="store_true", dest="suggest",
            help="Suggest reasonable parameters for aperture photometry.",
            default=False)

    return parser

def check_common_args(args,function_name,
                      valid_functions=['gaperture','gmap','gfind']):
    """Checks validity of some command line arguments used in gFind,
    gAperture, gMap, etc.  Returns the appropriate arguments as variables back
    to the calling procedure.
    """
    try:
        function_name = function_name.strip().lower()
    except AttributeError:
        raise gPhotonArgsError("Invalid function: {f}".format(f=function_name))

    if not function_name in valid_functions:
        raise gPhotonArgsError("Invalid function: {f}".format(f=function_name))

    try:
        args.band = args.band.strip().upper()
    except AttributeError:
        raise gPhotonArgsError("Invalid band: {b}".format(b=args.band))

    if not (args.ra or args.dec) and not args.skypos:
        raise gPhotonArgsError("Must specify either RA/DEC or SKYPOS.")
    elif (args.ra and args.dec) and args.skypos:
        raise gPhotonArgsError("Must specify either RA/DEC or SKYPOS, not both.")
    elif (args.ra and args.dec) and not args.skypos:
        args.skypos = [args.ra,args.dec]
    elif not (args.ra and args.dec) and args.skypos:
        args.ra, args.dec = args.skypos

    if args.suggest and function_name in ['gfind','gaperture']:
        (args.ra, args.dec, args.radius, args.annulus1,
            args.annulus2) = dbt.suggest_parameters(args.band, args.skypos,
                                                    retries=args.retries,
                                                    verbose=0)
        args.skypos = [args.ra, args.dec]
        if args.verbose:
            print "Recentering on ["+str(ra)+", "+str(dec)+"]"
            print "Setting radius to "+str(radius)
            print "Setting annulus to ["+str(annulus1)+", "+str(annulus2)+"]"

    if args.skypos:
        if np.array(args.skypos).shape != (2,):
            raise gPhotonArgsError(
                "Skypos (--skypos) must be a 2-element array.")
        args.ra, args.dec = args.skypos

    if args.ra and not (0. <= args.ra <= 360.):
        raise gPhotonArgsError(
            "RA of {ra} does not satisfy 0 <= RA <= 360".format(ra=args.ra))

    if args.dec and not (-90 <= args.dec <= 90):
        raise gPhotonArgsError(
            "Dec of {dec} does not satisfy -90 <= DEC <= 90".format(
                                                                dec=args.dec))

    if args.detsize and args.detsize <= 0.:
        raise gPhotonArgsError(
            "Effective field diameter (--detsize) must be > 0")

    if args.maxgap and args.maxgap <= 0.:
        raise gPhotonArgsError(
            "Maximum gap length (--maxgap) must be > 0 seconds.")

    if args.minexp and args.minexp <= 0.:
        raise gPhotonArgsError(
            "Minimum valid exposure depth (--minexp) must be > 0 seconds.")

    if args.retries and args.retries <= 0.:
        raise gPhotonArgsError("Number of retries (--retries) must be > 0.")

    # FIXME: This should accept a list of ranges...
    if args.trange:
        if len(args.trange) != 2:
            raise gPhotonArgsError(
                "trange (--trange) must be a 2-element list.")
        elif args.trange[0] <= 0. or args.trange[1] <= 0.:
            raise gPhotonArgsError(
                "trange (--trange) elements must both be > 0.")
        elif args.trange[0] >= args.trange[1]:
            raise gPhotonArgsError(
                "Maximum time must be greater than minimum time "+
                "in trange(--trange)")
    elif args.tmin and args.tmax:
        args.trange = [args.tmin,args.tmax]

    if args.tmin and args.tmin <= 0.:
        raise gPhotonArgsError("T0 (--t0) must be > 0.")

    if args.tmax and args.tmax <= 0.:
        raise gPhotonArgsError("T1 (--t1) must be > 0.")

    if args.tmin >= args.tmax:
        raise gPhotonArgsError(
            "Minimum time (--t0) must be < maximum time (--t1).")

    return args
