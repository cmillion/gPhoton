# Relative response:
#  The "canonical" GALEX pipeline generated response maps by _stamping_ the flat
#  field image onto a sky map at 1 second intervals, recentered and rolled
#  according to the aspect solution (boresight pointing) at that time and scaled
#  for the effective exposure time over that 1 second interval. (Note that the
#  aspect solution has "steps" of 1 second.) This "low resolution" relative
#  response map was then scaled 4x into a "high resolution" relative response
#  map (rrhr) with the same pixel size as the count map (cnt) from which it was
#  divided to produce the intensity map (int, or fully calibrated image).
#  This method is reproduced in rrhr() in imagetools.py.
#
#  The canonical method for creation of the rrhr includes a number of
#  interpolations which are both computationally slow and potentially introduce
#  unaccountable error. To speed up creation of calibrated light curves,
#  gAperture uses a method where it finds the location on the flat that
#  corresponds to the location on the detector where the source of interest
#  falls (by inverting the aspect correction). It then averages over an
#  aperture matching what is being used for photometry, scales it by the
#  effective exposure time, and uses that as the response. As with the canonical
#  method, this happens once per second corresonding to the resolution of the
#  aspect correction.
#  This method is reproduced in aperture_response() in curvetools.py.

