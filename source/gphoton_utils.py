# read and plot functionality for gPhoton .csv lightcurve files as created by gAperture.
import os
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def read_lc(ifile):
    """Reads a light curve CSV into a pandas DataFrame"""
    ## Read in the data as a pandas DataFrame object.
    try:
        data_frame = pd.io.parsers.read_csv(ifile, skipinitialspace=True, names=["t0", "t1", "ap_radius", "exptime", "cps", "cpserr", "flux", "flux_err", "mag", "mag_err", "r_inner", "r_outer", "bkg", "response", "counts", "apcorr1", "apcorr2"])
        ## Calculate the timestamp in JD.
        data_frame["JD"] = ((data_frame["t1"] - data_frame["t0"]) / 2. + data_frame["t0"]) / 86400. + 2440587.5
        
        ## Replace the default index with the timestamp in JD.
        data_frame.set_index("JD", drop=True, inplace=True, verify_integrity=True)

        ## Return the pandas objects.
        return data_frame
    except IOError:
        print "*** Error:  Could not find the file " + ifile + "."

def plot_lc(data_frame):
    """Plots a lightcurve from a CSV file
    data_frame - pandas DataFrame from read_lc()
    """
    plt.plot(data_frame.index.values, data_frame["flux"], "ko")
    plt.show()

