#!/usr/bin/python

import ast
import argparse
from regtestutils import datamaker
import MCUtils as mc
import gFind
import numpy as np
import os

def find_random_positions(rarange=[0,360],decrange=[-90,90],nsamples=10,
                          seed=323):
    np.random.seed(seed=seed%2)
    ra = np.random.uniform(rarange[0],rarange[1],nsamples)
    np.random.seed(seed=seed%3)
    dec = np.random.uniform(decrange[0],decrange[1],nsamples)
    return ra, dec

def calrun(outfile,band,nsamples=10,seed=323,rarange=[0,360],decrange=[-90,90],
           exprange=[0.,5000.],maglimit=24,verbose=0):
    """Generate a bunch of magnitudes with comparisons against MCAT values for
    random points on the sky within given legal ranges. Write it to a CSV.
    """
    (ra, dec) = find_random_positions(rarange=rarange,decrange=decrange,
                                      nsamples=nsamples,seed=seed)
    if verbose:
        print 'Running {n} random samples with seed of {seed}.'.format(
            n=nsamples,seed=seed)
        print 'Bounded by RA:[{r0},{r1}] and Dec:[{d0},{d1}]'.format(
            r0=rarange[0],r1=rarange[1],d0=decrange[0],d1=decrange[1])
        print 'Actual positions used will be:'
        print '{pos}'.format(pos=zip(ra,dec))

    for skypos in zip(ra,dec):
        expt = gFind.gFind(skypos=skypos,band=band,quiet=True)[band]['expt']
        if exprange[0]<=expt<=exprange[1]:
            print skypos, expt, True
            datamaker(band,skypos,outfile,maglimit=maglimit)
        else:
            print skypos, expt, False
    return

def setup_parser():
    parser = argparse.ArgumentParser(description="Generate a bunch of "+
        "magnitudes with comparisons against MCAT values for random points on "+
        "the sky within given legal ranges. Write it to a CSV.")
    parser.add_argument("-f", "--file", action="store", type=str, dest="file",
		default=None, help="File name (full path) for the output CSV.",
        required=True)
    parser.add_argument("-b", "--band", action="store", type=str.upper,
        dest="band", help="Band of NUV or FUV.", default='FUV', required=True,
        choices=["NUV","FUV", "nuv", "fuv"])
    parser.add_argument("-n", "--nsamples", action="store", type=int,
        dest="nsamples", help="Number of random locations to draw from sky.",
        default=10)
    parser.add_argument("--seed", action="store", type=int, dest="seed",
        help="Seed for the random number generator -- for reproducibility.",
        default=323)
    parser.add_argument("--rarange", action="store", dest="rarange",
		type=ast.literal_eval, default=[0.,360.],
        help="Two element list denoting valid ra range.")
    parser.add_argument("--decrange", action="store", dest="decrange",
		type=ast.literal_eval, default=[-90.,90.],
        help="Two element list denoting valid dec range.")
    parser.add_argument("--exprange", action="store", dest="exprange",
    	type=ast.literal_eval, default=[0.,5000.],
        help="Two element list denoting valid observation depths (in seconds).")
    parser.add_argument("--maglimit", action="store", type=float, default=24.,
        dest="maglimit", help="Lower limit of MCAT magnitudes to use.")
    parser.add_argument("-v", "--verbose", action="store", type=int, default=0,
        dest="verbose", help="Level of verbosity.", choices=[0,1,2])

    return parser

def check_args(args):
    if not (0<=args.rarange[0]<=360 and 0<=args.rarange[1]<=360 and
        args.rarange[0]<args.rarange[1]) or not (-90<=args.decrange[0]<=90 and
        -90<=args.decrange[1]<=90 and args.decrange[0]<args.decrange[1]):
        raise SystemExit("Invalid RA or Dec ranges.")
    if not (args.exprange[0]<args.exprange[1]):
        raise SystemExit("Invalid exposure range: {r}".format(r=args.exprange))
    return args

def __main__():
    args = setup_parser().parse_args()
    args = check_args(args)
    calrun(args.file,args.band,nsamples=args.nsamples,seed=args.seed,
           rarange=args.rarange,decrange=args.decrange,exprange=args.exprange,
           maglimit=args.maglimit,verbose=args.verbose)

if __name__ == "__main__":
    try:
        __main__()
    except (KeyboardInterrupt, pycurl.error):
        exit('Received Ctrl + C... Exiting! Bye.', 1)
