"""
.. module:: gphoton_utils

   :synopsis: Read and plot functionality for gPhoton .csv lightcurve files as
   created by gAperture.

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import galextools as gt

# ------------------------------------------------------------------------------
def read_lc(csvfile, comment='|'):
    """
    Read a light curve csv file from gAperture.

    :param csvfile: The name of the csv file to read.

    :type csvfile: str

    :param comment: The character used to denote a comment row.

    :type comment: str

    :returns: pandas DataFrame -- The contents of the csv file.
    """

    return pd.io.parsers.read_csv(csvfile, comment=comment)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def plot_lc(data_frame):
    """
    Plots a lightcurve from a CSV file data_frame - pandas DataFrame from
    read_lc()
	"""

    plt.plot(data_frame.index.values, data_frame["flux"], "ko")

    plt.show()

    # @CHASE - Don't need a return here?@
    return
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def model_errors(catmag, band, sigma=3., mode='mag', trange=[1, 1600]):
    """
    Give upper and lower expected bounds as a function of the nominal
	magnitude of a source. Very useful for identifying outliers.

    :param catmag: Nominal AB magnitude of the source.

    :type catmag: float

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param sigma: How many sigma out to set the bounds.

    :type sigma: float

    :param mode: Units in which to report bounds. Either 'cps' or 'mag'.

    :type mode: str

    :param trange: Set of integration times to compute the bounds on, in
    seconds.

    :type trange: list

    :returns: tuple -- A two-element tuple containing the lower and upper
    bounds, respectively.
	"""

    if mode != 'cps' and mode != 'mag':
        print 'mode must be set to "cps" or "mag"'
        exit(0)

    x = np.arange(trange[0], trange[1])

    cnt = gt.mag2counts(catmag, band)

    ymin = (cnt*x/x)-sigma*np.sqrt(cnt*x)/x

    ymax = (cnt*x/x)+sigma*np.sqrt(cnt*x)/x

    if mode == 'mag':
        ymin = gt.counts2mag(ymin, band)
        ymax = gt.counts2mag(ymax, band)

    return ymin, ymax
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def data_errors(catmag, t, band, sigma=3., mode='mag'):
    """
    Given an array (of counts or mags), return an array of 1-sigma error values.

    :param catmag: Nominal AB magnitude of the source.

    :type catmag: float

    :param t: Set of integration times to compute the bounds on, in seconds.

    :type t: list @CHASE - is this scalar or list? Also, consider trange
    instead of t to match first method?@

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param sigma: How many sigma out to set the bounds.

    :type sigma: float

    :param mode: Units in which to report bounds. Either 'cps' or 'mag'.

    :type mode: str

    :returns: tuple -- A two-element tuple containing the lower and upper
    uncertainty, respectively.
    """

    if mode != 'cps' and mode != 'mag':
        print 'mode must be set to "cps" or "mag"'
        exit(0)

    cnt = gt.mag2counts(catmag, band)

    ymin = (cnt*t/t)-sigma*np.sqrt(cnt*t)/t

    ymax = (cnt*t/t)+sigma*np.sqrt(cnt*t)/t

    if mode == 'mag':
        ymin = gt.counts2mag(ymin, band)
        ymax = gt.counts2mag(ymax, band)

    return ymin, ymax
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def dmag_errors(t, band, sigma=3., mode='mag', mags=np.arange(13, 24, 0.1)):
    """
    Given an exposure time, give dmag error bars at a range of magnitudes.

    :param t: Set of integration times to compute the bounds on, in seconds.

    :type t: list @CHASE - is this scalar or list? Also, consider trange
    instead of t to match first method?@

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param sigma: How many sigma out to set the bounds.

    :type sigma: float

    :param mode: Units in which to report bounds. Either 'cps' or 'mag'.

    :type mode: str

    :param mags: Set of magnitudes to compute uncertainties on.

    :type mags: numpy.ndarray

    :returns: tuple -- A three-element tuple containing the magnitudes and
    their lower and upper uncertainties, respectively.
    """

    cnts = gt.mag2counts(mags, band)*t

    ymin = (cnts-sigma/np.sqrt(cnts))/t

    ymax = (cnts+sigma/np.sqrt(cnts))/t

    if mode == 'mag':
        ymin = mags-gt.counts2mag(ymin, band)
        ymax = mags-gt.counts2mag(ymax, band)

    return mags, ymin, ymax
# ------------------------------------------------------------------------------
