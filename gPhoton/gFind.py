#!/usr/bin/python

"""
.. module:: gFind
   :synopsis: Find total amount of available data by position.
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
import argparse
# gPhoton imports.
import gPhoton.dbasetools as dbt
import gPhoton.gphoton_args as gargs

# ------------------------------------------------------------------------------
def gfind(band='both', detsize=1.1, exponly=False, gaper=False, maxgap=1500.0,
          minexp=1.0, quiet=False, retries=100, skypos=None, trange=None,
          verbose=0, skyrange=None):
    """
    Primary program in the module. Prints time ranges to the screen and
        returns the total exposure time as a float.

    :param band: The band being used, either 'BOTH', 'FUV', or 'NUV'.

    :type band: str

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :param exponly: If True, only report the exposure times.

    :type exponly: bool

    :param gaper: Return time ranges in a format that can be copy-pasted as
        a valid gAperture call.

    :type gaper: bool

    :param maxgap: Maximum gap size, in seconds, for data to be considered
        contiguous.

    :type maxgap: float

    :param minexp: Minimum gap size, in seconds, for data to be considered
        contiguous.

    :type minexp: float

    :param quiet: If True, don't print anything to STDOUT. Overrides verbose.

    :type quiet: bool

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param trange: Minimum and maximum time range to make a light curve,
        in GALEX time.

    :type trange: list

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param skyrange: RA and Dec extents, in degrees, defining the lengths of
        sides of a box on the sky that circumscribes the region of interest.

    :type skyrange: list

    :returns: dict -- The available data for the requested band(s).
	"""

    # Determine if we have to loop over both bands or just one.
    if band.upper() == 'BOTH':
        output = {'NUV':None, 'FUV':None}
    elif band.upper() in ['NUV', 'FUV']:
        output = {band.upper():None}
    else:
        raise SystemExit('Invalid band: {b}'.format(b=band))

    for this_band in list(output.keys()):
        # Get valid time ranges, but only if trange is not provided.
        ranges = dbt.fGetTimeRanges(this_band, skypos, maxgap=maxgap,
                                    minexp=minexp, verbose=verbose,
                                    detsize=detsize, trange=trange,
                                    skyrange=skyrange)
        if not ranges.any():
            if not quiet:
                print('No {band} exposure'
                      ' time in database.'.format(band=this_band))
            output[this_band] = {'expt':0, 't0':None, 't1':None}
        else:
            expt = (ranges[:, 1]-ranges[:, 0]).sum()
            if not quiet:
                print("{band}: {expt}s (raw) in {n} exposures.".format(
                    band=this_band, expt=expt, n=len(ranges)))
            if not exponly:
                if gaper:
                    f = '['
                    for r in ranges:
                        f += '[%.3f' % r[0] + ', %.3f' % r[1] + '],'
                    if not quiet:
                        print(f[:-1]+']')
                else:
                    for r in ranges:
                        if not quiet:
                            print('    [ %.3f' % r[0] + ', %.3f' % r[1] +
                                  ' ], %.3f' % (r[1]-r[0]) + ' seconds')
            output[this_band] = {'expt':expt, 't0':ranges[:, 0],
                                 't1':ranges[:, 1],
                                 'nearest_source':dbt.find_nearest_mcat(
                                     this_band, skypos, 0.05)}

    return output
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def setup_parser(iam='gfind', parser=None):
    """
    If a parser object is not provided, make one here.

    :param iam: The name of the method being called.

    :type iam: str

    :param parser: The command-line arguments.

    :type parser: argparse.ArgumentParser Namespace

    :returns: argparse.ArgumentParser Namespace -- The command-line arguments.
    """

    if parser is None:
        parser = argparse.ArgumentParser(description="Locate available data.")

    parser = gargs.common_args(parser, iam)
    parser.add_argument("--alt", "--gaper", action="store_true",
                        dest="gaper", default=False,
                        help="Format the output so that it can be copied and"
                        " pasted directly into a gMap or gAperture command"
                        " line?")
    parser.add_argument("--total", "--exponly", action="store_true",
                        dest="exponly", default=False, help="Report only the"
                        " total raw exposure time available in the database.")
    parser.add_argument("--quiet", action="store_true", dest="quiet",
                        help="Suppress all information to STDOUT.",
                        default=False)

    return parser
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def check_args(args, iam='gfind', allow_no_coords=False):
    """
    Checks validity of command line arguments and, in some cases,
        modifies them a little bit.

    :param args: The command-line arguments.

    :type args: argparse.ArgumentParser Namespace

    :param iam: The name of the method being called.

    :type iam: str

    :param allow_no_coords: If True, do not require coordinates to be provided.

    :type allow_no_coords: bool

    :returns: argparse.ArgumentParser Namespace -- The updated command-line
        arguments.
    """

    args = gargs.check_common_args(args, iam, allow_no_coords=allow_no_coords)

    return args
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def __main__():
    """
    Called when gFind is executed directly through the command line.
    """

    args = setup_parser().parse_args()
    args = check_args(args)

    exp_times = gfind(band=args.band, detsize=args.detsize,
                      exponly=args.exponly, gaper=args.gaper,
                      maxgap=args.maxgap, minexp=args.minexp, quiet=args.quiet,
                      retries=args.retries, skypos=args.skypos,
                      trange=args.trange, verbose=args.verbose,
                      skyrange=args.skyrange)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
if __name__ == "__main__":
    try:
        __main__()
    except KeyboardInterrupt:
        exit('Received Ctrl + C... Exiting! Bye.', 1)
# ------------------------------------------------------------------------------
