"""
.. module:: galextools
   :synopsis: This module contains helper and reference functions that are
       specific to GALEX and used by other modules, but do not directly access
       the database or generate end products of any kind.

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
from astropy import wcs as pywcs
import datetime
import time
import numpy as np

# ------------------------------------------------------------------------------
GPSSECS = 315532800+432000
# ------------------------------------------------------------------------------

def recovery_tranges():
    """
    Defines and returns an array of time ranges during which the spacecraft
        was in some sort of recovery mode (e.g. FUV cycling or CSP) and
        therefore any data from these periods should be viewed skeptically
        (because of things like observing while not at HVNOM).
    """
    return [# CSP (05-04-2010 to 06-23-2010)
        [time.mktime(datetime.date(2010, 5, 4).timetuple())-GPSSECS,
         time.mktime(datetime.date(2010, 6, 23).timetuple())-GPSSECS],
        ]

# ------------------------------------------------------------------------------
def isPostCSP(t, switch=961986575.):
    """
    Given a GALEX time stamp, return TRUE if it corresponds to a "post-CSP"
        eclipse. The actual CSP was on eclipse 37423, but the clock change
        (which matters more for calibration purposes) occured on 38268
        (t~=961986575.)

    :param t: The time stamp to test.

    :type t: float

    :param switch: The GALEX time stamp that defines pre- and post-CSP.

    :type switch: float

    :returns: bool -- Does this time correspond to a post-CSP eclipse?
    """

    # Check that the tscale has been applied properly
    if not switch/100. < t < switch*100.:
        raise ValueError('Did you apply tscale wrong?')

    return t >= switch
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def zpmag(band):
    """
    Define the zero point magnitude offset for the APER MCAT values.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :returns: float -- The zero point magnitude offset.
    """

    return {'NUV':20.08, 'FUV':18.82}[band]
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def aper2deg(apercode):
    """
    Convert SExtractor APER numbers to decimal degrees radii.

    :param apercode: The SExtractor APER number to convert.

    :type apercode: int

    :returns: float -- The APER radii in decimal degrees.
    """

    if not apercode == int(apercode) or apercode < 1 or apercode > 7:
        print("Error: `apercode` must be an integer in interval [1,7].")
        return None

    apers = np.array([1.5, 2.3, 3.8, 6.0, 9.0, 12.8, 17.3, 30., 60., 90.])/3600.

    return apers[int(apercode)-1]
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def apcorrect1(radius, band):
    """
    Compute an apeture correction. First way. Uses the table data in Figure 4
        from Morissey, et al., 2007

    :param radius: The photometric radius, in degrees.

    :type radius: float

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :returns: float -- The aperture correction.
    """

    # [Future]: Handle arrays.

    if not band in ['NUV', 'FUV']:
        print("Invalid band.")
        return

    aper = np.array([1.5, 2.3, 3.8, 6.0, 9.0, 12.8, 17.3, 30., 60., 90.])/3600.

    if radius > aper[-1]:
        return 0.

    if band == 'FUV':
        dmag = [1.65, 0.96, 0.36, 0.15, 0.1, 0.09, 0.07, 0.06, 0.03, 0.01]
    else:
        dmag = [2.09, 1.33, 0.59, 0.23, 0.13, 0.09, 0.07, 0.04, -0.00, -0.01]
        if radius > aper[-2]:
            return 0.

    if radius < aper[0]:
        return dmag[0]

    ix = np.where((aper-radius) >= 0.)
    x = [aper[ix[0][0]-1], aper[ix[0][0]]]
    y = [dmag[ix[0][0]-1], dmag[ix[0][0]]]
    m, C = np.polyfit(x, y, 1)

    return m*radius+C
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def apcorrect2(radius, band):
    """
    Compute an aperture correction in mag based upon an aperture radius in
        degrees. Second way. Uses the data in Table 1 from
        http://www.galex.caltech.edu/researcher/techdoc-ch5.html

    :param radius: The photometric radius, in degrees.

    :type radius: float

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :returns: float -- The aperture correction.
    """

    # [Future]: Handle arrays.

    if not band in ['NUV', 'FUV']:
        print("Invalid band.")
        return

    aper = np.array([1.5, 2.3, 3.8, 6.0, 9.0, 12.8, 17.3])/3600.

    if band == 'FUV':
        dmag = [1.65, 0.77, 0.2, 0.1, 0.07, 0.05, 0.04]
    else:
        dmag = [1.33, 0.62, 0.21, 0.12, 0.08, 0.06, 0.04]

    if radius > aper[-1]:
        # [Future]: Fix this, it isn't quite right...
        return dmag[-1]

    if radius < aper[0]:
        return dmag[0]

    ix = np.where((aper-radius) >= 0.)
    x = [aper[ix[0][0]-1], aper[ix[0][0]]]
    y = [dmag[ix[0][0]-1], dmag[ix[0][0]]]
    m, C = np.polyfit(x, y, 1)

    return m*radius+C
# ------------------------------------------------------------------------------

#
# Conversion facts for the following five functions can be found here:
# http://galexgi.gsfc.nasa.gov/docs/galex/FAQ/counts_background.html
#

# ------------------------------------------------------------------------------
def photometric_repeatability(cps, expt, band):
    """
    Estimate the photometric repeatability vs. magnitude .
        See: http://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html

    :param cps: Source flux in counts per second.

    :type cps: float

    :param expt: Effective exposure time, in seconds.

    :type expt: float

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :returns: float -- Estimated photometric repeatability based on flux.
    """

    scale = 0.050 if band == 'FUV' else 0.027

    return -2.5*(np.log10(cps)-
                 np.log10(cps+np.sqrt(cps*expt+(scale*cps*expt)**2.)/expt))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def detbg(area, band):
    """
    Nominal background in counts per second per 1.5" pixel.
        See: http://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html

    :param area: The area to calculate the background in square degrees.

    :type area: float

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :returns: float -- The nominal background for the given band, in counts per
        second.
    """

    rate = 1e-4 if band == 'FUV' else 1e-3

    return area * rate / ((1.5/(60*60))**2.)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def counts2mag(cps, band):
    """
    Converts GALEX counts per second to AB magnitudes.
        See: http://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html

    :param cps: The flux in counts per second.

    :type cps: float

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :returns: float -- The converted flux in AB magnitudes.
    """

    scale = 18.82 if band == 'FUV' else 20.08

    with np.errstate(invalid='ignore'):
        mag = -2.5 * np.log10(cps) + scale

    return mag
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def mag2counts(mag, band):
    """
    Converts AB magnitudes to GALEX counts per second.
        See: http://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html

    :param mag: The AB magnitude to convert.

    :type mag: float

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :returns: float -- The converted flux in counts per second.
    """

    scale = 18.82 if band == 'FUV' else 20.08

    return 10.**(-(mag-scale)/2.5)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def counts2flux(cps, band):
    """
    Converts GALEX counts per second to flux (erg sec^-1 cm^-2 A^-1).
        See: http://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html

    :param cps: The flux in counts per second.

    :type cps: float

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :returns: float -- The converted flux in erg sec^-1 cm^-2 A^-1.
    """

    scale = 1.4e-15 if band == 'FUV' else 2.06e-16

    return scale*cps
# ------------------------------------------------------------------------------

#
# END methods that use conversion factors from Goddard.
#

def local_nl_correction(mr, band):
    """
    Measured counts per second to predicted counts per sectond.
        Attempts to correct a measured count rate for nonlinearity per the
        formula given in Fig. 8 of Morrissey 2007.
    """

    coeffs = {'NUV':[-0.314, 1.365, -0.103],
              'FUV':[-0.531, 1.696, -0.225]}
    C0, C1, C2 = coeffs[band]

    return 10**np.roots([C2, C1, C0-np.log10(mr)])[1]

def deg2pix(skypos, skyrange, pixsz=0.000416666666666667):
    """
    Converts degrees to GALEX pixels rounded up to the nearest pixel
        so that the number of degrees specified will fully fit into the frame.

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param skyrange: Values in degrees RA and Dec of a box around skypos that
        defines the extent of the region of interest.

    :type skyrange: list

    :param pixsz: Width of a GALEX pixel, in degrees.

    :type pixsz: float

    :returns: float -- The converted number of pixels.
    """

    # [Future]: Fix this.
    # > gt.deg2pix([0,90],[0,0])
    # >  array([ 0.,  1.])
    # > gt.deg2pix([0,0],[0,0])
    # >  array([ 1.,  0.])
    wcs = pywcs.WCS(naxis=2)
    wcs.wcs.cdelt = [pixsz, pixsz]
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    # Set the reference pixel to [0,0] (which is actually the default)
    wcs.wcs.crpix = [0., 0.]

    # Fix the lower bound [RA, Dec] to the reference pixel
    wcs.wcs.crval = [skypos[0]-skyrange[0]/2., skypos[1]-skyrange[1]/2.]

    # Find the pixel position of the upper bound [RA,Dec]
    coo = [skypos[0]+skyrange[0]/2., skypos[1]+skyrange[1]/2.]

    # This is the image dimensions (at this location)
    # [Future]: Fix this: Is it better to ceil() or to floor() here?
    # Because ceil() completely captures the range but the actual
    # values of the edge pixels are suspect. Whereas floor() might
    # not catch the whole range, but all the values are meaningful.
    return np.abs(
        np.floor(wcs.sip_pix2foc(wcs.wcs_world2pix([coo], 1), 1)[0]))[::-1]
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def flat_scale_parameters(band):
    """
    Return the flat scaling parameters for a given band.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :returns: tuple -- A three-element tuple containing the flat scaling
        parameters.
    """

    if band == 'NUV':
        flat_correct = -0.0154
        flat_t0 = 840418779.02
        flat_correct_0 = 1.9946352
        flat_correct_1 = -1.9679445e-09
        flat_correct_2 = 9.3025231e-19
    elif band == 'FUV':
        flat_correct = -0.0031
        flat_t0 = 840418779.02
        flat_correct_0 = 1.2420282
        flat_correct_1 = -2.8843099e-10
        flat_correct_2 = 0.000
    else:
        print("Band not specified.")
        exit(1)

    # It turns out that flat_correct and flat_t0 never get used.
    # They are only retained above for historical purposes.
    return flat_correct_0, flat_correct_1, flat_correct_2
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def compute_flat_scale(t, band, verbose=0):
    """
    Return the flat scale factor for a given time.
        These are empirically determined linear scales to the flat field
        as a function of time due to diminished detector response. They
        were determined by Tom Barlow and are in the original GALEX pipeline
        but there is no published source of which I am aware.

    :param t: Time stamp(s) to retrieve the scale factor for.

    :type t: numpy.ndarray

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :returns: numpy.ndarray
    """

    if verbose:
        print("Calculating flat scale for t=", t, ", and band=", band)

    (flat_correct_0, flat_correct_1,
     flat_correct_2) = flat_scale_parameters(band)

    t = np.array(t)

    flat_scale = flat_correct_0+(flat_correct_1*t)+(flat_correct_2*t)*t

    # There's a bulk shift in the response after the CSP
    ix = np.where(t >= 881881215.995)

    if len(ix[0]):
        try:
            flat_scale[ix] *= 1.018
        except (TypeError, IndexError):
            # If it does not have '__getitem__' (because it's a scalar)
            flat_scale *= 1.018 if t >= 881881215.995 else 1.

    if verbose:
        print("         flat scale = ", flat_scale)

    return flat_scale
# ------------------------------------------------------------------------------
