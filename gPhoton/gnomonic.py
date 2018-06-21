"""
.. module:: gnomonic
   :synopsis: This module contains gnomonic projection methods for translating
       between detector and sky coordinates.
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
import numpy as np

# ------------------------------------------------------------------------------
R2D = 180./3.141592658979
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def gnomrev_simple(xi, eta, ra0, dec0, crota, cdelt, cenpix):
    """
    A reverse gnomonic projection.

    :param xi: Detector-space coordinte.

    :type xi: numpy.ndarray

    :param eta: Detector-space coordinate.

    :type eta: numpy.ndarray

    :param ra0: The right ascension of the boresite, in degrees.

    :type ra0: float

    :param dec0: The declination of the boresite, in degrees.

    :type dec0: float

    :param crota: The field rotation.

    :type crota: float

    :param cdelt: The resolution per pixel element at the field center.

    :type cdelt: float

    :param cenpix: The coordinates of the center pixel.

    :type cenpix: float

    :returns: tuple -- A two-element tuple containing the right ascension and
        declination, in degrees.
    """

    x = (xi - cenpix)*cdelt/R2D
    y = (eta - cenpix)*cdelt/R2D
    crotar = crota/R2D
    coscrotar = np.cos(crotar)
    sincrotar = np.sin(crotar)
    xp = x*coscrotar + y*sincrotar
    yp = -x*sincrotar + y*coscrotar
    rsq = xp*xp + yp*yp
    r = np.sqrt(rsq)

    cosdelta = np.ones(len(rsq))
    ix = (rsq > 1e-16).nonzero()
    cosdelta[ix] = 1. / np.sqrt(rsq[ix] + 1)

    sindelta = np.sqrt(1. - cosdelta*cosdelta)

    cosbeta, sinbeta = np.ones(len(r)), np.zeros(len(r))
    ix = (r > 0).nonzero()
    cosbeta[ix] = -yp[ix]/r[ix]
    sinbeta[ix] = -xp[ix]/r[ix]

    cosdec0 = np.cos(dec0/R2D)
    sindec0 = np.sin(dec0/R2D)
    xx = sindec0*sindelta*cosbeta + cosdec0*cosdelta
    yy = sindelta*sinbeta
    dec = np.arcsin(sindec0*cosdelta - cosdec0*sindelta*cosbeta)*R2D
    ra = ra0 + np.arctan2(yy, xx)*R2D
    ix = (ra >= 360).nonzero()
    ra[ix] -= 360
    ix = (ra < 0)
    ra[ix] += 360

    return ra, dec
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def gnomfwd_simple(ra, dec, ra0, dec0, crota, cdelt, cenpix):
    """
    A forward gnomonic projection.

    :param ra: The right ascension of the event, in degrees.

    :type ra: float numpy.ndarray

    :param dec: The declination of the event, in degrees.

    :type dec: numpy.ndarray

    :param ra0: The right ascension of the boresight.

    :type ra0: float

    :param dec0: The declination of the boresight.

    :type dec0: float

    :param crota: Rotation of the FOV.

    :type crota: float

    :param cdelt: Resolution at the center of the FOV.

    :type cdelt: float

    :param cenpix: Coordinates of the center pixel.

    :type cenpix: float

    :returns: tuple -- A two-element tuple containing the detector coordinates.
    """

    cosdec0 = np.cos(dec0/R2D)
    sindec0 = np.sin(dec0/R2D)
    drar = (ra0-ra)/R2D
    decr = dec/R2D
    cosdecr = np.cos(decr)
    sindecr = np.sin(decr)
    crotar = crota/R2D
    a = cosdecr * np.cos(drar)
    f = 1.0/(sindec0 * sindecr + a*cosdec0)
    xp = f * cosdecr * np.sin(drar)
    yp = f * (cosdec0 * sindecr - a*sindec0)
    cos_crotar = np.cos(-crotar)
    sin_crotar = np.sin(-crotar)
    x = xp*cos_crotar + yp*sin_crotar
    y = -xp*sin_crotar + yp*cos_crotar
    xi = x*R2D/cdelt + cenpix
    eta = y*R2D/cdelt + cenpix

    return xi, eta
# ------------------------------------------------------------------------------
