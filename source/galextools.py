# This module contains helper and reference functions that are specific
#  to GALEX and used by other modules but do not directly access the
#  database or generate end products of any kind.
import numpy as np
import math

gpssecs = 315532800+432000

# TODO: handle arrays
# Compute an aperture correction in mag based upon an aperture radius in
#  degrees. Uses the table data in Figure 4 from Morissey, et al., 2007
def apcorrect1(radius,band):
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

# Estimate the photometric repeatability vs. magnitude
def photometric_repeatability(cps,expt,band):
	scale = 0.050 if band=='FUV' else 0.027
	return -2.5 * (np.log10(cps) - np.log10(cps + np.sqrt(cps * expt + (scale * cps * expt)**2.)/expt))

# Computes the nominal detector background for a region
#  'area' is in square decimal degrees
def detbg(area, band):
	# Nominal background in counts per second per 1.5" pixel
	rate = 1e-4 if band=='FUV' else 1e-3
	return area * rate / ((1.5/(60*60))**2.)

# Converts GALEX counts per second to AB magnitudes
def counts2mag(cps,band):
	scale = 18.82 if band=='FUV' else 20.08
	return -2.5 * np.log10(cps) + scale

# Does the inverse of counts2mag
def mag2counts(mag,band):
	scale = 18.82 if band=='FUV' else 20.08
	return 10.**(-(mag-scale)/2.5)

# Converts GALEX counts per second to flux (erg sec^-1 cm^-2 A^-1)
def counts2flux(cps,band):
	scale = 1.4e-15 if band == 'FUV' else 2.06e-16
	return scale*cps

#
# END
#

# This converts degrees to GALEX pixels.
def deg2pix(degrees,CDELT2=0.000416666666666667):
        # Rounded up to the nearest pixel so that the number of
        #  degrees specified will fully fit into the frame.
        return math.ceil(degrees/CDELT2)

# These are empirically determined linear scales to the flat field
#  as a function of time due to diminished detector response. They
#  were determined by Tom Barlow and are in the original GALEX pipeline
#  but there is no published source of which I am aware.
def compute_flat_scale(t,band,verbose=1):
	if verbose:
		print "Calculating flat scale for t=",t,", and band=",band
	if band=='NUV':
		flat_correct=-0.0154
		flat_t0=840418779.02
		flat_correct_0=1.9946352
		flat_correct_1=-1.9679445e-09
		flat_correct_2=9.3025231e-19
	elif band=='FUV':
		flat_correct=-0.0154
		flat_t0=840418779.02
		flat_correct_0=1.2420282
		flat_correct_1=-2.8843099e-10
		flat_correct_2=0.000
	else:
		print "Band not specified."

	flat_scale=flat_correct_0+(flat_correct_1*t)+(flat_correct_2*t)*t

	if verbose:
		print "         flat scale = ",flat_scale

	return flat_scale

