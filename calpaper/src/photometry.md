##Photometry
The canonical GALEX pipeline used the Source Extractor (SExtractor) toolkit to extract photometric information from images. gAperture uses naive aperture photometry. These methods needs to be compared. The most straightforward way to do this is to compare gAperture flux measurements to their corresponding MCAT or GCAT values for a large number of sources across the sky. Do this efficiently will require enhancement of gAperture with the capability to determine and automatically use "optimum" aperture parameters based upon the MCAT/GCAT entry for the target in question.

A list of potential White Dwarf Standards appears in Table 5 of Morrissey, Patrick, et al. "The calibration and data products of GALEX." The Astrophysical Journal Supplement Series 173.2 (2007): 682.

###LDS749B
This is the primary white dwarf calibration standard for GALEX. It has nominal AB magnitudes of 14.71 in NUV and 15.57 in FUV.

`./gAperture.py --skypos [323.06766667,0.25400000] -a 0.005 --annulus [0.01,0.02] -b 'NUV' -v 2 -f 'LDS749B.csv' --maxgap 1600 --minexp 1000 --hrbg`


