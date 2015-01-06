# This module contains helper and reference functions that are specific
#  to GALEX and used by other modules but do not directly access the
#  database or generate end products of any kind.
import numpy as np
import math
from astropy import wcs as pywcs

GPSSECS = 315532800+432000

def zpmag(band):
    """Define the zero point magnitude offset for the APER MCAT values."""
    return {'NUV':20.08238,'FUV':18.81707}[band]

def aper2deg(aper):
    """Convert SExtractor APER numbers to decimal degrees radii."""
    if not aper==int(aper) or aper < 1 or aper > 7:
        print "APER must be an integer in interval [1,7]."
        return -1
    apers = np.array([1.5,2.3,3.8,6.0,9.0,12.8,17.3,30.,60.,90.])/3600.
    return apers[int(aper)-1]

# TODO: handle arrays
def apcorrect1(radius,band):
    """Compute an apeture correction. 1st way.
    radius - in degrees
    Uses the table data in Figure 4 from Morissey, et al., 2007
    """
    if not band in ['NUV','FUV']:
        print "Invalid band."
        return
    aper = np.array([1.5,2.3,3.8,6.0,9.0,12.8,17.3,30.,60.,90.])/3600.
    if radius>aper[-1]:
        return 0.
    if band=='FUV':
        dmag = [1.65,0.96,0.36,0.15,0.1,0.09,0.07,0.06,0.03,0.01]
    else:
        dmag = [2.09,1.33,0.59,0.23,0.13,0.09,0.07,0.04,-0.00,-0.01]
        if radius>aper[-2]:
            return 0.
    if radius<aper[0]:
        return dmag[0]

    ix  = np.where((aper-radius)>=0.)
    x   = [aper[ix[0][0]-1],aper[ix[0][0]]]
    y   = [dmag[ix[0][0]-1],dmag[ix[0][0]]]
    m,C = np.polyfit(x,y,1)
    return m*radius+C

# TODO: handle arrays
# Compute an aperture correction in mag based upon an aperture radius in
#  degrees. Uses the data in Table 1 from
#  www.galex.caltech.edu/research/techdoch-ch5.html
def apcorrect2(radius,band):
    """Compute an aperture correction. 2nd way.
    radius - in degrees
	Uses the data in Table 1 from
    www.galex.caltech.edu/research/techdoch-ch5.html
    """
    if not band in ['NUV','FUV']:
        print "Invalid band."
        return

    aper = np.array([1.5,2.3,3.8,6.0,9.0,12.8,17.3])/3600.
    if band=='FUV':
        dmag = [1.65,0.77,0.2,0.1,0.07,0.05,0.04]
    else:
        dmag = [1.33,0.62,0.21,0.12,0.08,0.06,0.04]
    if radius>aper[-1]:
        # FIXME: This isn't quite right...
        return dmag[-1]
    if radius<aper[0]:
        return dmag[0]

    ix  = np.where((aper-radius)>=0.)
    x   = [aper[ix[0][0]-1],aper[ix[0][0]]]
    y   = [dmag[ix[0][0]-1],dmag[ix[0][0]]]
    m,C = np.polyfit(x,y,1)
    return m*radius+C

#
# Conversion facts for the following five functions can be found here:
#  http://galexgi.gsfc.nasa.gov/docs/galex/FAQ/counts_background.html
#

def photometric_repeatability(cps,expt,band):
    """Estimate the photometric repeatability vs. magnitude."""
    scale = 0.050 if band=='FUV' else 0.027
    return -2.5*(np.log10(cps)-
                 np.log10(cps+np.sqrt(cps*expt+(scale*cps*expt)**2.)/expt))

def detbg(area, band):
	"""Nominal background in counts per second per 1.5" pixel"""
	rate = 1e-4 if band=='FUV' else 1e-3
	return area * rate / ((1.5/(60*60))**2.)

def counts2mag(cps,band):
	"""Converts GALEX counts per second to AB magnitudes"""
	scale = 18.82 if band=='FUV' else 20.08
	return -2.5 * np.log10(cps) + scale

def mag2counts(mag,band):
	"""Converts AB magnitudes to GALEX counts per second."""
	scale = 18.82 if band=='FUV' else 20.08
	return 10.**(-(mag-scale)/2.5)

def counts2flux(cps,band):
	"""Converts GALEX counts per second to flux (erg sec^-1 cm^-2 A^-1)"""
	scale = 1.4e-15 if band == 'FUV' else 2.06e-16
	return scale*cps

#
# END
#

def deg2pix(skypos,skyrange,pixsz=0.000416666666666667):
	"""Converts degrees to GALEX pixels rounded up to the nearest pixel
	so that the number of degrees specified will fully fit into the frame.
	"""
    # FIXME
    # > gt.deg2pix([0,90],[0,0])
    # >  array([ 0.,  1.])
    # > gt.deg2pix([0,0],[0,0])
    # >  array([ 1.,  0.])
	wcs = pywcs.WCS(naxis=2) # NAXIS = 2
	wcs.wcs.cdelt = [pixsz,pixsz]
	wcs.wcs.ctype = ['RA---TAN','DEC--TAN']
	# Set the reference pixel to [0,0] (which is actually the default)
	wcs.wcs.crpix = [0.,0.]
	# Fix the lower bound [RA, Dec] to the reference pixel
	wcs.wcs.crval = [skypos[0]-skyrange[0]/2.,skypos[1]-skyrange[1]/2.]
	# Find the pixel position of the upper bound [RA,Dec]
	coo = [skypos[0]+skyrange[0]/2.,skypos[1]+skyrange[1]/2.]
	# This is the image dimensions (at this location)
	# FIXME: Is it better to ceil() or to floor() here?
	#		 Because ceil() completely captures the range but the actual
	#		 values of the edge pixels are suspect. Whereas floor() might
	#		 not catch the whole range, but all the values are meaningful.
	return np.abs(np.floor(wcs.sip_pix2foc(wcs.wcs_world2pix([coo],1),1)[0]))[::-1]

def flat_scale_parameters(band):
    """Return the flat scaling parameters for a given band."""
    if band=='NUV':
        flat_correct = -0.0154
        flat_t0 = 840418779.02
        flat_correct_0 = 1.9946352
        flat_correct_1 = -1.9679445e-09
        flat_correct_2 = 9.3025231e-19
    elif band=='FUV':
        flat_correct = -0.0031
        flat_t0 = 840418779.02
        flat_correct_0 = 1.2420282
        flat_correct_1 = -2.8843099e-10
        flat_correct_2 = 0.000
    else:
        print "Band not specified."
        exit(1)
    # It turns out that flat_correct and flat_t0 never get used.
    # They are only retained above for historical purposes.
    return flat_correct_0, flat_correct_1, flat_correct_2

def compute_flat_scale(t,band,verbose=0):
    """Return the flat scale factor for a given time.
    These are empirically determined linear scales to the flat field
    as a function of time due to diminished detector response. They
    were determined by Tom Barlow and are in the original GALEX pipeline
    but there is no published source of which I am aware.
    """
    if verbose:
        print "Calculating flat scale for t=",t,", and band=",band
    (flat_correct_0, flat_correct_1,
                                flat_correct_2) = flat_scale_parameters(band)
    t = np.array(t)
    flat_scale=flat_correct_0+(flat_correct_1*t)+(flat_correct_2*t)*t
    # There's a bulk shift in the response after the CSP
    ix = np.where(t>=881881215.995)
    if len(ix[0]):
        try:
            flat_scale[ix] *= 1.018
        except TypeError:
            # if it does not have '__getitem__' (because it's a scalar)
            flat_scale *= 1.018 if t>=881881215.995 else 1.

    if verbose:
        print "         flat scale = ",flat_scale

    return flat_scale
