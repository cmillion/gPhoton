#!/usr/bin/python

"""
.. module:: gCalrun
   :synopsis: Performs a batch run of photometric extractions on semi-random
       sets of known sources for calibration / regression purposes.

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
from builtins import zip
import argparse
import ast
import numpy as np
# gPhoton imports.
import gPhoton.galextools as gt
from gPhoton import gFind
from gPhoton.regtestutils import datamaker

# ------------------------------------------------------------------------------
def find_random_positions(rarange=[0., 360.], decrange=[-90., 90.], nsamples=10,
                          seed=323):
    """
    Generate 'nsamples' random positions uniformly selected from the specified
        RA and DEC ranges.

    :param rarange: The minimum and maximum RA range to sample from.

    :type rarange: list

    :param decrange: The minimum and maximum DEC range to sample from.

    :type decrange: list

    :param nsamples: The number of random positions to return.

    :type nsamples: int

    :param seed: The seed to use when generating the random sample.

    :type seed: int

    :returns: tuple -- A two-element tuple containing 'nsamples' of RA and DEC
        values.
    """

    np.random.seed(seed=seed%2)

    ra = np.random.uniform(rarange[0], rarange[1], nsamples)

    np.random.seed(seed=seed%3)

    dec = np.random.uniform(decrange[0], decrange[1], nsamples)

    return ra, dec
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def calrun(outfile, band, nsamples=10, seed=323, rarange=[0., 360.],
           decrange=[-90., 90.], exprange=[0., 5000.], maglimit=24., verbose=0,
           radius=gt.aper2deg(4), annulus=[0.0083, 0.025]):
    """
    Generate a bunch of magnitudes with comparisons against MCAT values for
        random points on the sky within given legal ranges. Write it to a CSV.

    :param outfile: Name of the output file.

    :type outfile: str

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param nsamples: Number of random positions to sample.

    :type nsamples: int

    :param seed: The seed to use when generating the random sample.

    :type seed: int

    :param rarange: The minimum and maximum RA range to sample from.

    :type rarange: list

    :param decrange: The minimum and maximum DEC range to sample from.

    :type decrange: list

    :param exprange: The minimum and maximum exposure time to sample,
        in seconds.

    :type exprange: list

    :param maglimit: The faintest source to consider.

    :type maglimit: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param radius: Photometric aperture radius, in degrees.

    :type radius: float

    :param annulus: Inner and outer extent of background annulus, in degrees.

    :type annulus: list
    """

    (ra, dec) = find_random_positions(rarange=rarange, decrange=decrange,
                                      nsamples=nsamples, seed=seed)

    if verbose:
        print('Running {n} random samples with seed of {seed}.'.format(
            n=nsamples, seed=seed))
        print('Bounded by RA:[{r0},{r1}] and Dec:[{d0},{d1}]'.format(
            r0=rarange[0], r1=rarange[1], d0=decrange[0], d1=decrange[1]))
        print('Actual positions used will be:')
        print('{pos}'.format(pos=list(zip(ra, dec))))

    for skypos in zip(ra, dec):
        expt = gFind(skypos=skypos, band=band, quiet=True)[band]['expt']
        if exprange[0] <= expt <= exprange[1]:
            print(skypos, expt, True)
            datamaker(band, skypos, outfile, maglimit=maglimit, verbose=verbose,
                      searchradius=0.01)
        else:
            print(skypos, expt, False)

    return
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
def setup_parser():
    """
    Setup input arguments.
    """

    parser = argparse.ArgumentParser(
        description="Generate a bunch of "
        "magnitudes with comparisons against MCAT values for random points on "
        "the sky within given legal ranges. Write it to a CSV.")

    parser.add_argument("-a", "--aperture", action="store", dest="radius",
                        help="Aperture radius in decimal degrees", type=float)
    parser.add_argument("-i", "--inner", action="store", dest="annulus1",
                        help="Inner annulus radius for background subtraction",
                        type=float)
    parser.add_argument("-o", "--outer", action="store", dest="annulus2",
                        help="Outer annulus radius for background subtraction",
                        type=float)
    parser.add_argument("--annulus", action="store", type=ast.literal_eval,
                        dest="annulus", help="Annulus inner and outer radius"
                        " definition as '[inner,outer]'")
    parser.add_argument("-f", "--file", action="store", type=str, dest="file",
                        default=None, help="File name (full path) for the"
                        " output CSV.", required=True)
    parser.add_argument("-b", "--band", action="store", type=str.upper,
                        dest="band", help="Band of NUV or FUV.", default='FUV',
                        required=True, choices=["NUV", "FUV", "nuv", "fuv"])
    parser.add_argument("-n", "--nsamples", action="store", type=int,
                        dest="nsamples", help="Number of random locations to"
                        " draw from sky.", default=10)
    parser.add_argument("--seed", action="store", type=int, dest="seed",
                        help="Seed for the random number generator -- for"
                        " reproducibility.", default=323)
    parser.add_argument("--rarange", action="store", dest="rarange",
                        type=ast.literal_eval, default=[0., 360.],
                        help="Two element list denoting valid ra range.")
    parser.add_argument("--decrange", action="store", dest="decrange",
                        type=ast.literal_eval, default=[-90., 90.],
                        help="Two element list denoting valid dec range.")
    parser.add_argument("--exprange", action="store", dest="exprange",
                        type=ast.literal_eval, default=[0., 5000.],
                        help="Two element list denoting valid observation"
                        " depths (in seconds).")
    parser.add_argument("--maglimit", action="store", type=float, default=24.,
                        dest="maglimit", help="Lower limit of MCAT magnitudes"
                        " to use.")
    parser.add_argument("-v", "--verbose", action="store", type=int, default=0,
                        dest="verbose", help="Level of verbosity.",
                        choices=[0, 1, 2])

    return parser
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def check_args(args):
    """
    Ensures that the input RA is between 0. and 360., and the DEC is between
        -90 and 90.

    :param args: The command-line arguments.

    :type args: argparse.ArgumentParser Namespace

    :returns: argparse.ArgumentParser Namespace -- The command-line arguments.
        Included for function usage consistency.
    """

    if (not(0. <= args.rarange[0] <= 360. and 0. <= args.rarange[1] <= 360. and
            args.rarange[0] < args.rarange[1]) or not (
                -90. <= args.decrange[0] <= 90. and -90. <= args.decrange[1] <=
                90. and args.decrange[0] < args.decrange[1])):
        raise SystemExit("Invalid RA or Dec ranges.")

    if not args.exprange[0] < args.exprange[1]:
        raise SystemExit("Invalid exposure range: {r}".format(r=args.exprange))

    return args
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def __main__():
    """
    Called when gCalrun is executed through the command line.
    """

    args = setup_parser().parse_args()
    args = check_args(args)

    calrun(args.file, args.band, nsamples=args.nsamples, seed=args.seed,
           rarange=args.rarange, decrange=args.decrange, exprange=args.exprange,
           maglimit=args.maglimit, verbose=args.verbose, radius=args.radius,
           annulus=args.annulus)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
if __name__ == "__main__":
    try:
        __main__()
    except (KeyboardInterrupt):
        exit('Received Ctrl + C... Exiting! Bye.', 1)
# ------------------------------------------------------------------------------
