"""
.. module:: gphoton_utils
   :synopsis: Read, plot, time conversion, and other functionality useful when
       dealing with gPhoton data.

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# gPhoton imports.
import gPhoton.galextools as gt

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
        print('mode must be set to "cps" or "mag"')
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
def data_errors(catmag, band, t, sigma=3., mode='mag'):
    """
    Given an array (of counts or mags), return an array of n-sigma error values.

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
        print('mode must be set to "cps" or "mag"')
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

# ------------------------------------------------------------------------------
def calculate_jd(galex_time):
    """
    Calculates the Julian date, in the TDB time standard, given a GALEX time.

    :param galex_time: A GALEX timestamp.

    :type galex_time: float

    :returns: float -- The time converted to a Julian date, in the TDB
    time standard.
    """

    if np.isfinite(galex_time):
        # Convert the GALEX timestamp to a Unix timestamp.
        this_unix_time = Time(galex_time + 315964800., format="unix",
                              scale="utc")

        # Convert the Unix timestamp to a Julian date, measured in the
        # TDB standard.
        this_jd_time = this_unix_time.tdb.jd
    else:
        this_jd_time = np.nan

    return this_jd_time
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def calculate_jd_utc(galex_time):
    """
    Calculates the Julian date, in the UTC time standard, given a GALEX time.

    :param galex_time: A GALEX timestamp.

    :type galex_time: float

    :returns: float -- The time converted to a Julian date, in the UTC
    time standard.
    """

    if np.isfinite(galex_time):
        # Convert the GALEX timestamp to a Unix timestamp.
        this_unix_time = Time(galex_time + 315964800., format="unix",
                              scale="utc")

        # Convert the Unix timestamp to a Julian date, measured in the
        # UTC standard.
        this_jd_time = this_unix_time.utc.jd
    else:
        this_jd_time = np.nan

    return this_jd_time
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def calculate_jd_tai(galex_time):
    """
    Calculates the Julian date, in the TAI time standard, given a GALEX time.

    :param galex_time: A GALEX timestamp.

    :type galex_time: float

    :returns: float -- The time converted to a Julian date, in the TAI
    time standard.
    """

    if np.isfinite(galex_time):
        # Convert the GALEX timestamp to a Unix timestamp.
        this_unix_time = Time(galex_time + 315964800., format="unix",
                              scale="utc")

        # Convert the Unix timestamp to a Julian date, measured in the
        # UTC standard.
        this_jd_time = this_unix_time.tai.jd
    else:
        this_jd_time = np.nan

    return this_jd_time
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def calculate_caldat(galex_time):
    """
    Calculates a Gregorian calendar date given a GALEX time, in the UTC
    time standard.

    :param galex_time: A GALEX timestamp.

    :type galex_time: float

    :returns: float -- The time converted to a Gregorian calendar date,
    in the UTC time standard.
    """

    if np.isfinite(galex_time):
        # Convert the GALEX timestamp to a Unix timestamp.
        this_unix_time = Time(galex_time + 315964800., format="unix",
                              scale="utc")

        # Convert the Unix timestamp to a calendar date, measured in the
        # UTC standard.
        this_caldat_time = this_unix_time.iso
    else:
        this_caldat_time = 'NaN'

    return this_caldat_time
# ------------------------------------------------------------------------------

