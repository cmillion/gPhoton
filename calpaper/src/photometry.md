##Photometry
The canonical GALEX pipeline used the Source Extractor (SExtractor) toolkit to extract photometric information from images. gAperture uses naive aperture photometry. These methods needs to be compared. The most straightforward way to do this is to compare gAperture flux measurements to their corresponding MCAT or GCAT values for a large number of sources across the sky. Do this efficiently will require enhancement of gAperture with the capability to determine and automatically use "optimum" aperture parameters based upon the MCAT/GCAT entry for the target in question.

A list of potential White Dwarf Standards appears in Table 5 of Morrissey, Patrick, et al. "The calibration and data products of GALEX." The Astrophysical Journal Supplement Series 173.2 (2007): 682.

Morrissey, et al. suggests a 12" (0.003 degree) aperture as a good generic value.

"For both detectors, the onset of saturation is at approximately mAB ~ 15 (corresponding to 34 counts s^-1 FUV and 108 counts s^-1 NUV) _in the lowest-gain regions of the detector_."

###LDS749B
This is the primary white dwarf calibration standard for GALEX. It has nominal AB magnitudes of *14.71* in NUV and *15.57* in FUV. Note that this is past the saturation point in NUV, which needs to be accounted for in the analysis.

`./gAperture.py --skypos [323.06766667,0.25400000] -a 0.005 --annulus [0.01,0.02] -b 'NUV' -v 2 -f 'LDS749B_NUV_rr.csv' --maxgap 1600 --minexp 1000 --hrbg --response --overwrite`

In [222]: counts2mag(mag2counts(data['mag']-data['apcorr2'],'NUV').mean(),'NUV')Out[222]: 14.713075354913997

`./gAperture.py --skypos [323.06766667,0.25400000] -a 0.005 --annulus [0.01,0.02] -b 'FUV' -v 2 -f 'LDS749B_FUV_rr.csv' --maxgap 1600 --minexp 1000 --hrbg --response --overwrite`

In [240]: counts2mag(mag2counts(data['mag']-data['apcorr2'],'FUV').mean(),'FUV')
Out[240]: 15.414742366262452


