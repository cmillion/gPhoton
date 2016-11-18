#!/usr/bin/python

"""
.. module:: gPipeline
   :synopsis: Runs the module that runs the standalone gPhoton calibration
       pipelines to create aspect-corrected, time-tagged photon events from
       low-level archived GALEX data products.

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
import argparse
from builtins import str
# gPhoton imports.
from gPhoton.PhotonPipe import photonpipe

# ------------------------------------------------------------------------------
def gpipeline(raw6file, scstfile, band, outbase, aspfile, ssdfile, nullout,
              retries=100):
    """
    Wrapper that calls the PhotonPipe method.

    :param raw6file: Name of the raw6 file to use.

    :type raw6file: str

    :param scstfile: Spacecraft state file to use.

    :type scstfile: str

    :param band: Name of the band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param outbase: Base of the output file names.

    :type outbase: str

    :param aspfile: Name of aspect file to use.

    :type aspfile: str

    :param ssdfile: Name of Stim Separation Data file to use.

    :type ssdfile: str

    :param nullfile: Name of output file to record NULL lines.

    :type nullfile: str

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int
    """

    photonpipe(raw6file, scstfile, band, outbase, aspfile, ssdfile, nullout,
               retries=retries)

    return
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def setup_parser():
    """
    Defines command-line arguments.

    :returns: argparse.ArgumentParser Namespace -- The command-line arguments.
    """

    parser = argparse.ArgumentParser(description="Generate time-tagged lists"
                                     " of aspect-corrected photons positions.")
    parser.add_argument("-r", "--raw6", action="store", type=str,
                        dest="raw6file", help="raw6 path/filename")
    parser.add_argument("-a", "--aspect", action="store", type=str,
                        dest="aspfile", help="aspection path/filename")
    parser.add_argument("-s", "--scst", action="store", type=str,
                        dest="scstfile", help="spacecraft state path/filename")
    parser.add_argument("-b", "--band", action="store", type=str,
                        dest="band", help="[NF]UV band designation")
    parser.add_argument("-o", "--outbase", action="store", type=str,
                        dest="outbase", help="output file(s) path/filename"
                        " base", default=None)
    parser.add_argument("-d", "--ssd", action="store", type=str,
                        dest="ssdfile", help="stim separation (SSD)"
                        " path/filename")
    parser.add_argument("-u", "--nullout", action="store_true", dest="nullout",
                        help="write NULL entries to a separate file",
                        default=False)
    parser.add_argument("--retries", action="store", type=int, default=20,
                        help="Query attempts before timeout.")

    return parser
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def check_args(args):
    """
    Checks validity of command line arguments and, in some cases
        modifies them a little bit.

    :param args: The command-line arguments.

    :type args: argparse.ArgumentParser Namespace

    :returns: argparse.ArgumentParser Namespace -- The updated command-line
        arguments.
    """

    if not args.raw6file:
        print("Must specify a RAW6 filename (--raw6).")
        raise SystemExit

    if not args.scstfile:
        print("Must specify a SCST filename (--scst).")
        raise SystemExit

    if not args.outbase:
        print("Must specify an output base filename (--outbase).")
        raise SystemExit

    if args.aspfile:
        args.aspfile = str(args.aspfile).split(',')

    if args.ssdfile:
        args.ssdfile = str(args.ssdfile)

	# If the band is not explicity called, attempt to derive it from the raw6
    # filename.
    if not args.band:
        print("Determining band from raw6 filename...")
        if '-fd-raw6' in args.raw6file:
            args.band = 'FUV'
        elif '-nd-raw6' in args.raw6file:
            args.band = 'NUV'
        else:
            print("Unable to parse band from raw6 filename. Specify band on"
                  " command line using --band.")
            raise SystemExit
    else:
        args.band = args.band.upper()
        if not args.band in ["NUV", "FUV"]:
            print("Band must be NUV or FUV. ")
            raise SystemExit

    return args
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def __main__():
    """
    Called when gPipeline is executed directly through the command line.
    """

    args = setup_parser().parse_args()

    args = check_args(args)

    gpipeline(args.raw6file, args.scstfile, args.band, args.outbase,
              args.aspfile, args.ssdfile, args.nullout, retries=args.retries)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
if __name__ == "__main__":
    try:
        __main__()
    except KeyboardInterrupt:
        exit('Received Ctrl + C... Exiting! Bye.', 1)
# ------------------------------------------------------------------------------
