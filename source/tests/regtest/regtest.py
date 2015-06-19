"""This is a full system regression test for gPhoton. It takes a while to run,
perhaps many hours. It will print progress as it goes and then a final verdict
of PASS/FAIL. You should run it from the terminal command line as a script like
> python regtest.py
"""
import os
import sys
import pandas as pd
sys.path=sys.path+['../..'] # Probably an awful hack that you should never do.
from gCalrun import calrun

print 'GENERATING TEST CSV DATA (may take a while)'
for band in ['NUV','FUV']:
    print 'BAND: {b}'.format(b=band)
    calrun('DB10_calrun_{b}_test.csv'.format(b=band),band,nsamples=1,seed=323,
        rarange=[0,360],decrange=[53,90],verbose=1,calpath='../../../cal/')

print 'BEGINNING REGRESSION TEST'
import more_test

print 'DELETING TEST CSV DATA'
for band in ['NUV','FUV']:
    os.system('rm DB10_calrun_{b}_test.csv'.format(b=band))
