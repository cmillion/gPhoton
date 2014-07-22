""" This module contains functions for checking arguments across gFind,
gAperture, gMap, and other gPhoton functions.
"""

import argparse
import os
import ast

class gPhotonArgsError(Exception):
    """Exception class specific to gphoton_args."""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def setup_args(function_name=None):
	"""Defines the arguments and options for the parser object when
	called from the command line.  Accepts a string used to determine which
    arguments to add to the Parser object.  Valid function names are "gfind",
    "gaperture", or "gmap" (all case insensitive).
	"""
        if function_name is not None:
            function_name = function_name.strip().lower()
            if (function_name != 'gfind' and function_name != 'gaperture'
                                         and function_name != 'gmap'):
                raise gPhotonArgsError("Choice of function not understood. Could not setup command-line arguments.  Name provided = "+str(function_name)+".")
        else:
            raise gPhotonArgsError("Name of function must be provided to setup command-line arguments.")

        # Initiate argument parser object with help description appropriate
        # for the function.
        if function_name == 'gfind':
            parser = argparse.ArgumentParser(description="Locate available data at specified coordinates and time intervals.")
        elif function_name == 'gaperture':
            parser = argparse.ArgumentParser(description="Create lightcurve CSV files at specified coordinates and time intervals.")
        elif function_name == 'gmap':
            parser = argparse.ArgumentParser(description="Create intensity maps and animated movie cubes at specified coordinates and time intervals.")
        else:
            # Should never get here, but raise an Exception just in case.
            raise gPhotonArgsError("Name of function not understood.  Name provided = " + function_name+".")

        # Add arguments to the Parser object per function.

        #ADDHDR is accepted by gAperture only.
        if function_name == 'gaperture':
            parser.add_argument("--addhdr", action="store_true", dest="addhdr",
            default=False,
            help="Write the command line and column headers to the top of the .csv file?  Default = False.")

        # ALT is accepted by gFind only.
        if function_name == 'gfind':
            parser.add_argument("--alt", "--gaper", action="store_true",
            dest="gaper", default=False, help="Format the output so that it can be copied and pasted directly into a gMap or gAperture command line? Default = False.")

        #ANGLE is accepted by gMap only.
        if function_name == 'gmap':
            parser.add_argument("--angle", action="store", type=float,
            dest="angle", default=None,
            help="The angle subtended in both RA and DEC in degrees.")

        # ANNULUS is accepted by gAperture only.
        if function_name == 'gaperture':
            parser.add_argument("--annulus", action="store", dest="annulus",
            help="Annulus inner and outer radius definition with format [inner,outer] in degrees.")

        # APERTURE is accepted by gAperture only.
        if function_name == 'gaperture':
            parser.add_argument("-a", "--aperture", action="store", type=float,
            dest="radius", help="Aperture radius in decimal degrees.")

        # Band is accepted by gFind, gAperture, and gMap, but only gFind can accept "both" as a valid choice.
        if function_name == 'gfind':
            parser.add_argument("-b", "--band", action="store", type=str,
            dest="band", help="Band designation, choice of {FUV, NUV, BOTH}.  Default is BOTH.",
            metavar="BAND", default="BOTH", 
            choices=["NUV","FUV","BOTH", "nuv", "fuv", "both"])
        elif function_name == 'gaperture' or function_name == 'gmap':
            parser.add_argument("-b", "--band", action="store", type=str,
            dest="band", help="Band designation, choice of {FUV, NUV}.",
            metavar="BAND", choices=["NUV","FUV","nuv","fuv"])

        #BESTPARAMS is accepted by gAperture only.
        if function_name == 'gaperture':
            parser.add_argument("--bestparams", "--best", action="store_true",
            dest="best", default=False,
            help="Set parameters to produce the highest quality lightcurve?  Potentially slow.  Default = False.")

        #CALPATH is accepted by gAperture and gMap.
        if function_name == 'gaperture' or function_name == 'gmap':
            parser.add_argument("--calpath", action="store", type=str,
            dest="calpath", default=os.pardir+os.sep+"cal"+os.sep,
            help="Path to the directory that contains the calibration files.  Default = '../cal/'")

        #COADD is accepted by gAperture and gMap only.
        if function_name == 'gaperture' or function_name == 'gmap':
            parser.add_argument("--coadd", action="store_true", dest="coadd",
            default=False,
            help="Return the coadded flux (gAperture) or a coadded image (gMap) over all requested time ranges?  Default = False.")

        #COUNT is accepted by gMap only.
        if function_name == 'gmap':
            parser.add_argument("--count", action="store", type=str,
            dest="cntfile", default=None, help="File name (full path) for the count image.")

        # DEC is accepted by gFind, gAperture, and gMap.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            parser.add_argument("-d", "--dec", action="store", type=float,
            dest="dec", metavar="DEC",
            help="Center Declination position in decimal degrees.  Must be 0 < DEC < 90.")

        #DECANGLE is accepted by gMap only.
        if function_name == 'gmap':
            parser.add_argument("--decangle", action="store", type=float,
            dest="decangle", default=None,
            help="The angle of sky in degrees that the declination subtends.  Overrides --angle.")

        # DETSIZE is accepted by gFind, gAperture, and gMap.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            parser.add_argument("--detsize", action="store", type=float,
            dest="detsize", default=1.25,
            help="Set the effective field diameter in degrees for the exposure search.  Default = 1.25.")

        # EXPONLY is accepted by gFind only.
        if function_name == 'gfind':
            parser.add_argument("--total", "--exponly", action="store_true",
            dest="exponly", default=False,
            help="Report only the total raw exposure time available in the database, not the individual time ranges? Default = False.")

        # FILE is accepted by gAperture only.
        if function_name == 'gaperture':
            parser.add_argument("-f", "--file", action="store", type=str,
            dest="file",
            help="File name (full path) for CSV lightcurve file.",default=None)

        #HRBG is accepted by gAperture only.
        if function_name == 'gaperture':
            parser.add_argument("--usehrbg", "--hrbg", action="store_true",
            dest="usehrbg", default=False, 
            help="Use the higher quality 'swiss cheese' background estimation method?  Default = False.")

        # INNER is accepted by gAperture only.
        if function_name == 'gaperture':
            parser.add_argument("-i", "--inner", action="store", type=float,
            dest="annulus1",
            help="Inner annulus radius for background subtraction in degrees.")

        #INTENSITY is accepted by gMap only.
        if function_name == 'gmap':
            parser.add_argument("--intensity", action="store", type=str,
            dest="intfile", default=None,
            help="File name (full path) for the intensity image.")

        #MASKRAD is accepted by gMap only.
        if function_name == 'gmap':
            parser.add_argument("--maskrad", action="store", type=float,
            dest="maskrad", default=1.,
            help="The radius in degrees at which detector events will be masked out.  Default = 1.")

	# MAXGAP is accepted by gFind, gAperture, and gMap.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            parser.add_argument("-g", "--gap", "--maxgap", action="store",
            type=float, dest="gap", default=1500.,
            help="Maximum gap size in seconds for data to be considered contiguous.  Default = 1500.")

        #MEMLIGHT is accepted by gMap only.
        if function_name == 'gmap':
            parser.add_argument("--memlight", action="store", type=float,
            dest="memlight", default=100.,
            help="Reduce server-side memory usage by requesting data in chunks of no more than this depth in seconds.  Default = 100.")

        # MINEXP is accepted by gFind, gAperture, and gMap.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            parser.add_argument("--minexp", action="store", type=float, dest="minexp", help="Minimum contiguous exposure in seconds for data to be reported.  Default = 1.", default=1.)
                
        # OUTER is accepted by gAperture only.
        if function_name == 'gaperture':
            parser.add_argument("-o", "--outer", action="store", type=float, dest="annulus2", help="Outer annulus radius for background subtraction in degrees.")

        #OVERWRITE is accepted by gAperture and gMap only.
        if function_name == 'gaperture' or function_name == 'gmap':
            parser.add_argument("--overwrite", "--ow", "--clobber", action="store_true", dest="overwrite", help="Overwrite existing output files?  Default = False.", default=False)

        # QUIET is accepted by gFind only.
        if function_name == 'gfind':
            parser.add_argument("--quiet", action="store_true", dest="quiet", help="Suppress all information to STDOUT? Default = False.", default=False)

        # RA is accepted by gFind, gAperture, and gMap.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            parser.add_argument("-r", "--ra", action="store", type=float, dest="ra", help="Center Right Ascension position in decimal degrees.  Must be 0 < RA < 360.", metavar="RA")

        #RAANGLE is accepted by gMap only.
        if function_name == 'gmap':
            parser.add_argument("--raangle", action="store", type=float, dest="raangle", help="The angle of sky in degrees that the right ascension subtends.  Overrides --angle.",default=None)

        #RESPONSE is accepted by gMap only.
        if function_name == 'gmap':
            parser.add_argument("--response", action="store", type=str, dest="rrfile", help="File name (full path) for the response image.", default=None)

        # RETRIES is accepted by gFind, gAperture, and gMap.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            parser.add_argument("--retries", action="store", type=int, dest="retries", help="Set the number of times to ping the server for a response before defining a query failure.  Default is 20, set to a large number if you expect, or want to allow, the query to take a long time without issuing a failure.", default=20)
                
	# SKYPOS is accepted by gFind, gAperture, and gMap.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            parser.add_argument("--skypos", action="store", dest="skypos", help="Alternate method for specifying sky position with format '[RA,Dec]'",type=ast.literal_eval)

        #STAMP is accepted by gAperture only.
        if function_name == 'gaperture':
            parser.add_argument("--stamp", action="store", type=str, dest="stamp", help="Output file name for JPEG preview image stamp of the targeted region.",default=None)

        # STEP is accepted by gAperture and gMap only.
        if function_name == 'gaperture' or function_name == 'gmap':
            parser.add_argument("-s", "--step", "--frame", action="store", type=float, dest="stepsz", help="Step size for lightcurve or movie in seconds.  Default = 0. (no binning).",default=0.)
        
        # SUGGEST is accepted by gFind and gAperture only.
        if function_name == 'gfind' or function_name == 'gaperture':
            parser.add_argument("--suggest", action="store_true", dest="suggest", help="Suggest optimum parameters for aperture photometry? Default = False.", default=False)

	# T0 is accepted by gFind, gAperture, and gMap.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            parser.add_argument("--t0", "--tmin", action="store", type=float, dest="tmin", help="Minimum date of observation to consider (specify in GALEX time standard).  Default = 1.",default=1.)

	# T1 is accepted by gFind, gAperture, and gMap.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            parser.add_argument("--t1", "--tmax", action="store", type=float, dest="tmax", help="Maxium date of observation to consider (specify in GALEX time standard).  Default = 1000000000000.",default=1000000000000.)

	# TRANGE is accepted by gFind, gAperture, and gMap.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            parser.add_argument("--trange", "--tranges", action="store", dest="trange", help="Time range in which to limit the search, in the format '[t0,t1]' (specify in GALEX time standard).", type=ast.literal_eval)

        #USERR is accepted by gAperture only.
        if function_name == 'gaperture':
            parser.add_argument("--userr", "--response", action="store_true", dest="userr", help="Use the relative response correction?  Default = False.", default=False)

	# VERBOSE is accepted by gFind, gAperture, and gMap.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            parser.add_argument("-v", "--verbose", action="store", type=int, dest="verbose", help="Prints extra information to STDOUT (higher number = more output). Choices are {0,1,2,3}, default = 0.", metavar="VRB",default=0.,choices=[0,1,2,3])

        # Return the Parser object.
	return parser


def check_args(args, function_name=None):
	"""Checks validity of some command line arguments used in gFind, gAperture, gMap, etc.  Returns the appropriate arguments as variables back to the calling procedure."""

        if function_name is not None:
            function_name = function_name.strip().lower()
            if function_name != 'gfind' and function_name != 'gaperture' and function_name != 'gmap':
                raise gPhotonArgsError("Choice of function not understood.  Could not check command-line arguments.  Name provided = " + function_name+".")
        else:
            raise gPhotonArgsError("Name of function must be provided to check command-line arguments.")

        # Check the angle argument.
        if function_name == 'gmap':
            angle = args.angle
            if angle is not None and angle <= 0.:
                raise gPhotonArgsError("Angle must be > 0.")

        # Check the annulus argument.
        if function_name == 'gaperture':
            annulus = args.annulus
            if annulus is not None:
                if isinstance(annulus, basestring):
                    annulus = ast.literal_eval(annulus)
                if len(annulus) != 2:
                    raise gPhotonArgsError("Annulus must be a 2-element list.")
                elif annulus[0] <= 0. or annulus[1] <= 0.:
                    raise gPhotonArgsError("Annulus inner and outer radii must both be > 0.")
                elif annulus[0] >= annulus[1]:
                    raise gPhotonArgsError("Inner annulus must be < outer annulus.")

        # Check the aperture argument.
        if function_name == 'gaperture':
            radius = args.radius
            if radius is not None and radius <= 0.:
                raise gPhotonArgsError("Aperture must be > 0.")

        # Check the BAND argument.  Convert it to all uppercase even if provided in lowercase.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            # Make sure band is in all uppercase and has whitespace removed.
            band = args.band.strip().upper()

        # Check the BEST argument.  If it is set, then USEHRBG and USERR both must be set to True as well.
        if function_name == 'gaperture':
            best=args.best
            if best:
                usehrbg = True
                userr = True
            else:
                usehrbg = args.usehrbg
                userr = args.userr

        # Check the COUNT argument.
        if function_name == 'gmap':
            cntfile = args.cntfile
            if cntfile is not None:
                file_base = os.path.dirname(cntfile)
                if file_base and not os.path.isdir(file_base):
                    os.makedirs(os.path.abspath(file_base,0755))            

        # Check the DEC argument.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            dec = args.dec
            if dec is not None and (dec < -90. or dec > 90.):
                raise gPhotonArgsError("DEC must be -90. <= DEC <= 90.")

        # Check the DECANGLE argument.
        if function_name == 'gmap':
            decangle = args.decangle
            if decangle is not None and decangle <= 0.:
                raise gPhotonArgsError("DEC angle must be > 0.")

        # Check the DETSIZE argument.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            detsize = args.detsize
            if detsize is not None and detsize <= 0.:
                raise gPhotonArgsError("Effective field diameter must be > 0.")

        # Check the FILE argument.
        if function_name == 'gaperture':
            outfile = args.file
            if outfile is not None:
                file_base = os.path.dirname(outfile)
                if file_base and not os.path.isdir(file_base):
                    os.makedirs(os.path.abspath(file_base,0755))            

        # Check the INNER argument.
        if function_name == 'gaperture':
            annulus1 = args.annulus1
            if annulus1 is not None and annulus1 <= 0.:
                raise gPhotonArgsError("Inner annulus must be > 0.")

        # Check the INTENSITY argument.
        if function_name == 'gmap':
            intfile = args.intfile
            if intfile is not None:
                file_base = os.path.dirname(intfile)
                if file_base and not os.path.isdir(file_base):
                    os.makedirs(os.path.abspath(file_base,0755))

        # Check the MASKRAD argument.
        if function_name == 'gmap':
            maskrad = args.maskrad
            if maskrad is not None and maskrad <= 0.:
                raise gPhotonArgsError("Mask radius must be > 0.")

        # Check the MAXGAP argument.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            maxgap = args.gap
            if maxgap is not None and maxgap <= 0.:
                raise gPhotonArgsError("Maximum gap length must be > 0 seconds.")

        # Check the MEMLIGHT argument.
        if function_name == 'gmap':
            memlight = args.memlight
            if memlight is not None and memlight <= 0.:
                raise gPhotonArgsError("Maximum data chunk per request must be > 0 seconds.")

        # Check the MINEXP argument.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            minexp = args.minexp
            if minexp is not None and minexp <= 0.:
                raise gPhotonArgsError("Minimum contiguous data size to report back must be > 0 seconds.")

        # Check the OUTER argument.
        if function_name == 'gaperture':
            annulus2 = args.annulus2
            if annulus2 is not None and annulus2 <= 0.:
                raise gPhotonArgsError("Outer annulus must be > 0.")

        # Check the RA argument.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            ra = args.ra
            if ra is not None and (ra < 0. or ra > 360.):
                raise gPhotonArgsError("RA must be 0. <= RA <= 360.")

        # Check the RAANGLE argument.
        if function_name == 'gmap':
            raangle = args.raangle
            if raangle is not None and raangle <= 0.:
                raise gPhotonArgsError("DEC angle must be > 0.")

        # Check the RESPONSE argument.
        if function_name == 'gmap':
            rrfile = args.rrfile
            if rrfile is not None:
                file_base = os.path.dirname(rrfile)
                if file_base and not os.path.isdir(file_base):
                    os.makedirs(os.path.abspath(file_base,0755))

        # Check the RETRIES argument.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            retries = args.retries
            if retries is not None and retries <= 0.:
                raise gPhotonArgsError("Number of retries must be > 0.")

        # Check the SKYPOS argument.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            skypos = args.skypos
            if skypos is not None:
                if isinstance(skypos, basestring):
                    skypos = ast.literal_eval(skypos)
                if len(skypos) != 2:
                    raise gPhotonArgsError("SKYPOS must be a 2-element list.")
                elif skypos[0] < 0. or skypos[0] > 360.:
                    raise gPhotonArgsError("RA must be 0. <= RA <= 360.")
                elif skypos[1] < -90. or skypos[1] > 90.:
                    raise gPhotonArgsError("DEC must be -90. <= DEC <= 90.")

        # Check the STAMP argument.
        if function_name == 'gaperture':
            stamp = args.stamp
            if stamp is not None:
                file_base = os.path.dirname(stamp)
                if file_base and not os.path.isdir(file_base):
                    os.makedirs(os.path.abspath(file_base,0755))

        # Check the STEP argument.
        if function_name == 'gaperture' or function_name == 'gmap':
            stepsz = args.stepsz
            if stepsz is not None and stepsz < 0.:
                raise gPhotonArgsError("Step size must be >= 0.")

        # Check the T0 argument.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            tmin = args.tmin
            if tmin is not None and tmin <= 0.:
                raise gPhotonArgsError("T0 (TMIN) must be > 0.")

        # Check the T1 argument.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            tmax = args.tmax
            if tmax is not None and tmax <= 0.:
                raise gPhotonArgsError("T1 (TMAX) must be > 0.")

        # Check the TRANGE argument.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            trange = args.trange
            if trange is not None:
                if isinstance(trange, basestring):
                    trange = ast.literal_eval(trange)
                if len(trange) != 2:
                    raise gPhotonArgsError("TRANGE must be a 2-element list.")
                elif trange[0] <= 0. or trange[1] <= 0.:
                    raise gPhotonArgsError("Time ranges must both be > 0.")
                elif trange[0] >= trange[1]:
                    raise gPhotonArgsError("Minimum time must be < maximum time in TRANGE.")
                
        #
        # Perform additional checks and actions below here.
        #

        # Make sure the user did not specify both STEPSZ and COADD
        if function_name == 'gaperture' or function_name == 'gmap':
            if stepsz and args.coadd:
                raise gPhotonArgsError("Cannot specify both STEPSZ and COADD.")

        # Make sure either the radius or suggest parameter is specified for gAperture.
        if function_name == 'gaperture':
            if not radius and not args.suggest:
                raise gPhotonArgsError("Must specify an aperture radius, or ask for a suggested radius.")

	# Print suggested parameters to screen if requested by gFind, or use them if requested by gAperture.
        if function_name == 'gfind' or function_name == 'gaperture':
            if args.suggest and function_name == 'gfind':
                out = suggest_parameters(band,skypos,verbose=1,retries=retries)
                print ""
            elif args.suggest and function_name == 'gaperture':
                ra,dec,radius,annulus1,annulus2 = suggest_parameters(band,skypos,retries=retries)
                if verbose:
                    print "Recentering on ["+str(ra)+", "+str(dec)+"]"
                    print "Setting radius to "+str(radius)
                    print "Setting annulus to ["+str(annulus1)+", "+str(annulus2)+"]"

	# Create the "skypos" list if skypos was not directly supplied on input.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            if not (ra and dec) and not skypos:
		raise gPhotonArgsError("Must specify either RA/DEC or SKYPOS.")
            elif (ra and dec) and skypos:
		raise gPhotonArgsError("Must specify either RA/DEC or SKYPOS, not both.")
            elif (ra and dec):
		skypos = [ra,dec]

        # Create the "skyrange" list if skyrange ("angle") was not directly supplied on input.
        if function_name == 'gmap':
            if angle is not None:
                skyrange = [angle,angle]
                if raangle is not None:
                    skyrange[0] = raangle
                if decangle is not None:
                    skyrange[1] = decangle
            elif raangle is not None and decangle is not None:
                skyrange = [raangle,decangle]
            else:
                raise gPhotonArgsError("Must specify either RAANGLE/DECANGLE or ANGLE.")

	# Make sure tmin < tmax, create the "trange" list if trange was not directly supplied on input, but tmin and tmax were.  Note that since there are default values we don't need to check if these variables are None.
        if function_name == 'gfind' or function_name == 'gaperture' or function_name == 'gmap':
            if tmin >= tmax:
                raise gPhotonArgsError("Minimum time (TMIN) must be < maximum time (TMAX).")
            if not trange:
		trange = [tmin,tmax]

        # Make sure the inner annulus is < outer annulus, create the "annulus" list if annulus was not directly supplied on input.
        if function_name == 'gaperture':
            if (annulus1 and annulus2) and annulus1 >= annulus2:
                raise gPhotonArgsError("Inner annulus must be < outer annulus.")
            if not (annulus1 and annulus2) and not annulus:
		raise gPhotonArgsError("Must specify either INNER/OUTER annuli or ANNULUS.")
            elif (annulus1 and annulus2) and annulus:
		raise gPhotonArgsError("Must specify either INNER/OUTER annuli or ANNULUS, not both.")
            elif (annulus1 and annulus2):
		annulus = [annulus1,annulus2]

        # If the user has selected to use HRBG, then annulus must be defined.
        if function_name == 'gaperture':
            if usehrbg and not annulus:
                raise gPhotonArgsError("You must specify a background annulus with HRBG.")

        # If the OVERWRITE option is not set, and output files already exist, then exit early to not waste computation on stuff that won't get saved.  We can consider changing this later to ask the user if they want to continue anyways, but I don't think that makes a lot of sense.
        if function_name == 'gaperture':
            if outfile:
                if not args.overwrite and os.path.exists(outfile):
                    raise gPhotonArgsError("Output lightcurve file "+outfile+" already exists.  Use --overwrite or specify a different output file name.")
            if stamp:
                if not args.overwrite and os.path.exists(stamp):
                    raise gPhotonArgsError("Output preview image file "+stamp+" already exists.  Use --overwrite or specify a different output file name.")
            
        if function_name == 'gmap':
            if cntfile:
                if not args.overwrite and os.path.exists(cntfile):
                    raise gPhotonArgsError("Count image file "+cntfile+" already exists.  Use --overwrite or specify a different output file name.")
            if intfile:
                if not args.overwrite and os.path.exists(intfile):
                    raise gPhotonArgsError("Intensity image file "+intfile+" already exists.  Use --overwrite or specify a different output file name.")
            if rrfile:
                if not args.overwrite and os.path.exists(rrfile):
                    raise gPhotonArgsError("Response image file "+rrfile+" already exists.  Use --overwrite or specify a different output file name.")
        
        # Return the appropriate set of parameters given the input function.  We do not return parameters that aren't used later, or that have no possibility of changing from their input values (e.g., Booleans that do not depend on other parameters).
	if function_name == 'gfind':
            return (band, detsize, maxgap, minexp, retries, skypos, trange)
        elif function_name == 'gaperture':
            return (annulus,band,best,detsize,outfile,maxgap,minexp,radius,retries,skypos,stamp,stepsz,trange,usehrbg,userr)
        elif function_name == 'gmap':
            return (band,cntfile,detsize,intfile,rrfile,skypos,maxgap,memlight,minexp,retries,skyrange,stepsz,trange)
