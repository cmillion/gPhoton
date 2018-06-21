#!/usr/bin/python

"""
.. module:: gMap
   :synopsis: Create images and image cubes.
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
import argparse
import os
import numpy as np
# gPhoton imports.
import gPhoton.dbasetools as dbt
import gPhoton.gphoton_args as gargs
from gPhoton.imagetools import write_images

# ------------------------------------------------------------------------------
def gmap(band, cntfile=None, coadd=None, detsize=1.1, intfile=None,
         skypos=None, maxgap=1500., memlight=100., minexp=1.,
         overwrite=False, retries=100, skyrange=None, stepsz=0., trange=None,
         verbose=0, cntcoaddfile=False, intcoaddfile=False):
    """
    Use a mix of strings (if we want to make an output file) and Booleans
        (False if we do not). I don't think this is the best way (I'd rather
        have separate variables for the True/False to create an image, and a
        string with the name of the output image), but for now this is kept
        since this is how the code was originally written.

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param cntfile: Name of the count file to make.

    :type cntfile: str

    :param coadd: If true, create a coadd image using all available data within
        the specified time range (don't use time bins).

    :type coadd: bool

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :param intfile: Name of intensity file to make.

    :type intfile: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param maxgap: Maximum gap size, in seconds, for data to be considered
        contiguous.

    :type maxgap: float

    :param memlight: If specified, breaks the queries into sections of this
        many seconds.

    :type memlight: float

    :param minexp: Minimum gap size, in seconds, for data to be considered
        contiguous.

    :type minexp: float

    :param overwrite: If True, overwite an existing output file.

    :type overwrite: bool

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int

    :param skyrange: RA and Dec extents, in degrees, that define the lengths
        of edges of a box on the sky, centered at skypos, that defines the
        sky region of interest.

    :type skyrange: list

    :param stepsz: The size of the time bins to use, in seconds.

    :type stepsz: float

    :param trange: Minimum and maximum time range to make images for,
        in GALEX time.

    :type trange: list

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param cntcoaddfile: Name of the count coadd file to make.

    :type cntcoaddfile: str

    :param intcoaddfile: Name of the intensity coadd file to make.

    :type intcoaddfile: str
    """

    # [Future]: Consider improving this section.
    write_cnt = cntfile if (cntfile) else False
    write_int = intfile if (intfile) else False
    write_cnt_coadd = cntcoaddfile if cntcoaddfile else False
    write_int_coadd = intcoaddfile if intcoaddfile else False

    # If gMap is called via an import, and no trange is specified,
    # then use a database call to get the appropriate time range.
    if trange is None:
        trange = dbt.fGetTimeRanges(band, skypos, trange=[6.E8, 11.E8],
                                    maxgap=maxgap, minexp=minexp,
                                    detsize=detsize, skyrange=skyrange)
    if np.array(trange).size == 0:
        print('No data available.')
        return

    if len(np.array(trange).shape) == 1:
        trange = [trange]

    # If trange is out of time order, this will put it back into order.
    trange = np.array(trange)[np.argsort(np.array(trange)[:, 0])]

    write_images(band.upper(), skypos, trange, skyrange, write_cnt=write_cnt,
                 write_int=write_int, framesz=stepsz,
                 overwrite=overwrite, verbose=verbose, memlight=memlight,
                 coadd=coadd, retries=retries, write_cnt_coadd=write_cnt_coadd,
                 write_int_coadd=write_int_coadd, detsize=detsize)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def setup_parser(iam='gmap'):
    """
    Defines command line arguments.

    :param iam: The name of the method being called.

    :type iam: str

    :returns: argparse.ArgumentParser Namespace -- The command-line arguments.
    """

    parser = argparse.ArgumentParser(description="Generate images / maps.")
    parser = gargs.common_args(parser, iam)
    parser.add_argument("--angle", "-a", action="store", type=float,
                        dest="angle", default=None,
                        help="The angle subtended in both RA and DEC in"
                        " degrees.")
    parser.add_argument("--raangle", action="store", type=float, dest="raangle",
                        help="The angle of sky in degrees that the right"
                        " ascension subtends. Overrides --angle.", default=None)
    parser.add_argument("--decangle", action="store", type=float,
                        dest="decangle", default=None,
                        help="The angle of sky in degrees that the declination"
                        " subtends. Overrides --angle.")
    parser.add_argument("--count", action="store", type=str, dest="cntfile",
                        default=None, help="File name (full path) for the"
                        " count image.")
    parser.add_argument("--intensity", action="store", type=str, dest="intfile",
                        default=None, help="File name (full path) for the"
                        " intensity image.")
    parser.add_argument("--count_coadd", action="store", type=str,
                        dest="cntcoaddfile", default=None,
                        help="File name (full path) for the count coadd image.")
    parser.add_argument("--intensity_coadd", action="store", type=str,
                        dest="intcoaddfile", default=None,
                        help="File name (full path) for the intensity coadd"
                        " image.")
    parser.add_argument("--memlight", action="store", type=float,
                        dest="memlight", default=100., help="Reduce server-"
                        "side memory usage by requesting data in chunks of"
                        " no more than this depth in seconds.")

    return parser
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def check_args(args, iam='gmap'):
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

    # Check that either skypos _or_ (raangle _and_ decangle) are defined.
    if ((not args.skyrange and not (args.raangle and args.decangle)
         and not args.angle)):
        raise SystemExit('Must specify either skyrange or ra/dec angle.')
    if ((args.raangle and not args.decangle) or (not args.raangle and
                                                 args.decangle)):
        raise SystemExit('Must specify both raangle and deangle or neither.')

    # --angle overwrites everything
    if args.angle:
        args.skyrange = [args.angle, args.angle]

    # --skryange overwrites everything else
    if args.skyrange:
        if np.array(args.skyrange).shape == (2,):
            args.raangle, args.decangle = args.skyrange
        else:
            gargs.gPhotonArgsError(
                "Invalid --skyrange: {s}".format(s=args.skyrange))

    # use --angle to fill in missing subtend angles
    if args.angle and not args.raangle:
        args.raangle = args.angle
    if args.angle and not args.decangle:
        args.decangle = args.angle

    # Check that requested image subtend angles are positive
    for a in [args.angle, args.raangle, args.decangle]:
        if not a:
            continue
        if a <= 0:
            raise SystemExit("Subtended angles must be > 0 degrees.")

    args.skyrange = [args.raangle, args.decangle]

    if args.detsize and args.detsize <= 0.:
        raise SystemExit("Mask radius must be > 0 degrees.")
    if args.memlight and args.memlight <= 0.:
        raise SystemExit("Maximum data chunk per must be > 0 seconds.")

    for image in [args.cntfile, args.intfile]:
        # Check for overwriting existing images.
        if image and os.path.exists(image) and not args.overwrite:
            raise SystemExit("{f} already exists.".format(f=image))

    return args
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def __main__():
    """
    Called when gMap is executed directly through the command line.
    """

    args = setup_parser().parse_args()

    args = check_args(args)

    if args.verbose:
        print('{a} image(s) at {pos}'.format(
            a='Writing' if not args.coadd else 'Coadding', pos=args.skypos))
        print('		of dimensions {w}x{h} degrees'.format(
            w=args.skyrange[0], h=args.skyrange[1]))
        print('		with a virtual detector {d} degrees across'.format(
            d=args.detsize))
        print('		in time range(s): {t}'.format(t=repr(args.trange)))

    gmap(band=args.band, cntfile=args.cntfile,
         coadd=args.coadd, detsize=args.detsize, intfile=args.intfile,
         skypos=args.skypos, maxgap=args.maxgap,
         memlight=args.memlight, minexp=args.minexp, overwrite=args.overwrite,
         retries=args.retries, skyrange=args.skyrange, stepsz=args.stepsz,
         trange=args.trange, verbose=args.verbose,
         cntcoaddfile=args.cntcoaddfile, intcoaddfile=args.intcoaddfile)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
if __name__ == "__main__":
    try:
        __main__()
    except KeyboardInterrupt:
        exit('Received Ctrl + C... Exiting! Bye.', 1)
# ------------------------------------------------------------------------------
