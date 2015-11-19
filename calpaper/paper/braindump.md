The GALEX spacecraft operated for some time over some years. The primary mission ended at a certain time followed by an extended mission, known as CAUSE, for some time after that.

The GALEX spacecraft observed in the near and far ultraviolet using two identically constructed microchannel plate detectors. Light from the boresight was redirected to the detectors by a dichroic and passed through filters.

The GALEX mission collected some number of observations, accounting for some amount of approximate total raw data volume and some other amount of approximate reduced data volume.

Microchannel plates are non-integrating photon detectors. That is, rather than accumulating detector events in integrated bins, the microchannel plates record position and time information for every event individually. Due to computer storage and processing constraints during the mission, the GALEX mission team only released integrated maps of these data on per-observation time scales. The high time and spatial resolution photon level data was never released to the community, except in isolated cases by special request (as _x-files_), and no serious attempts were made to understand the detector performance or calibration parameters for integrations shorter than about 200 seconds. Furthermore, by the end of the mission, the GALEX calibration pipeline software had grown to sufficient complexity that no attempts to get it to run outside of the GALEX network of computers at Caltech were successful. With the support of Space Telescope Science Institute, the authors undertook the project of building a standalone GALEX calibration pipeline--using it to archive the ~1.8 trillion GALEX photon events accumulated over the course of the whole mission--and a suite of software tools to make those photons usable at arbitrary spatial and temporal scales.

High time resolution astronomy is of increasing interest as observation and detector technology permit it. The GALEX data set represents a huge resource for this.

Engineering Considerations
1. Minimize external dependencies. Dependencies should be limited to Python packages considered to be mature or at least well established.
2. Eliminate mission bookkeeping conventions which arbitrarily separate data. A photon event is uniquely defined by time and position.
3. Design to the most likely use cases.

Calibration Challenges
1. The exposure time correction and, in particular, the global dead time correction was not reliable for short time ranges. This is because the dead time correction was computed by comparing the measured stime rate against the nominal stim rate of 79 counts per second. At short time ranges (on the order of seconds), the variance in this rate due simply to Poissonian counting statistics made the error in the exposure time calculation sufficiently high as to make flux calculations meaningless. The solution was to use the empirical dead time correction, a linear fitted relationship of dead time as a function of the global detector count rate (which can be thousands of counts per second). [from the calibration paper] Our tests have shown that this is a reliable estimator of the dead time correction and leads to stable/reliable effective exposure times at short time ranges.
2. The relative response (i.e. flat field) was not characterized at spatial resolutions that make it meaningful for short time domains. That is, because the detector is constantly in motion with respect to the sky during any observation, any given source samples a region of the detector. The resolution chosen for the flat field--about equal to the width of a single "dither"--was chosen with the assumption that an integrated observation would smear out any sub-resolution detector response effects, and the uncertainties would average out. This assumption does not hold for short exposure times. We have no solution for this; it continues to be a problem. But we know of efforts by other investigators to characterize the flat at higher spatial resolutions.

Uniqueness Properties of Photon Data
1. No need to integrate images. (No need to interpolate.)
2. Can directly sample the flat / mask. (No need to interpolate.)
3. Complete flexibility in defining time ranges and bins.
4. Data volume. There is a lot of it.
