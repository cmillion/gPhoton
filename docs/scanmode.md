##GALEX Scan Mode Data

###Overview
In the "post-NASA" or CAUSE phase of the GALEX mission, the spacecraft began observing frequently in a mode known as scan mode. In scan mode, the spacecraft boresight would traverse and observe many degrees of sky, rapidly, in a long swath. This was in contrast to the traditional boresight dither or the "petal pattern" and AIS modes in which the spacecraft hovered over single sky pointings and was reduced to a low voltage (non-observing) state when scanning orienting between pointings. The benefit of scan mode was that it allowed the team to more rapidly complete the All Sky Survey in the ultraviolet. But scan mode was performed at a time when the mission was operating on both minimal staff and budget, so while nominal data products were produced, the calibration pipeline was not optimized to handle this observing mode and no additional scientific support for the data was available or provided by the team.

###Considerations
The nature and peculiarities of scan mode data are still being teased out. Not very many people have looked closely at it, and no careful analysis has been done of how it differs from data in the traditional modes. This document will be updated as more is learned, and we welcome input from data users.

Scan mode has a high "FAIL" rate as set by both automated flagging in the calibration pipeline and manual QA. Some FAIL graded visits may contain useful information or be fine for specific investigative tasks. But in general, care and caution should be used when performing any investigation with scan mode data regardless of its PASS/FAIL status.

Below are some known or suspected _gotchas_ with respect to scan mode:
* The GALEX calibration pipeline assumed integrated images of a size and shape approximately the same as the detector FOV (a circle with a 1.25 degree diameter). To accomodate scan mode which would produce single observations covering many degrees, the calibration pipeline was hacked to subdivide each scan into multiple subvisits akin to AIS subvisits. This is the reason that the data is formatted this way in the products released by the mission.
* There is a high rate of FAIL vists. Some users have found that the PASS/FAIL status of a visit should be taken with a grain of salt. Some FAIL data appears to contain perfectly usable data and some PASS data has serious issues.
* A common artifact in the count and intensity maps is a doubling or ghosting of individual sources or images. We believe that this is caused by jumps in the refined aspect solution where the aspect correction stage of the calibration pipeline got confused.
* Photometry with GALEX _sort of_ assumes that the source of interest remains in the same general region of the detector over an observation (i.e. that the source is on the same part of the flat), and users are in general advised to not trust flux measurements near the edge of the detector. Both of these conditions are violated by scan mode where every source samples a cord across the whole detector which includes the edges. The method that the GALEX calibration pipeline used for trimming regions of questionable response--masking out portions of the image with an integrated response below a threshhold fraction of the maximum response--would not have worked at all for scan mode data in the downtrack direction.
* Coadds were not created for scan mode data.

###gPhoton
The GALEX photon database project at MAST, called gPhoton, has the potential to mitigate some of these known issues if not completely eliminate them. In particular, it should be possible to manually exclude regions with bad aspect solutions or data near the detector edge and generate coadded images and aperture photometry (regardless of mode or visit/sub-visit conventions). A project to improve the flat / response characterization or at least understand the photometric consequences of observing over wide swaths of the detector is in the works as well, though we welcome contributions to that effort.


