**User Documentation for the MAST/GALEX Photon Database and Tools**

Chase Million<sup>1</sup>, Bernie Shiao<sup>2</sup>, Scott Fleming<sup>2</sup>,
Myron Smith<sup>2</sup>

<sup>1</sup>Million Concepts (chase.million@gmail.com),
<sup>2</sup>Space Telescope Science Institute

###Summary

The MAST/GALEX photon database and tools exist in an effort to maximize the flexibility and utility of the GALEX data set. The GALEX detectors were microchannel plates which recorded detector position and time-of-arrival information for every detected photon event with a time resolution of five thousandths of a second, composing a huge and rich short time domain survey of the UV. Due to digital storage space and processing limitations, the data was only formally released by the mission team as integrated images. The "extended" photon list files--internally known as _x-files_--were only provided by special request and with little to no additional support for their calibration or use.

The official GALEX calibration pipeline software ("the canonical pipeline")--written in several languages with a sprawling network of dependencies--has also never been successfully ported to any system outside of the GALEX internal network at Caltech. Even though the source code for this pipeline will be made publicly available through the Mikulski Archive at Space Telescope (MAST) in the future, the end of the GALEX project would have effectively marked the end of the capability to generate photon level data specifically and revisit the GALEX calibration more generally. A suite of software tools called _gPhoton_ has been developed by the authors with support of MAST and Space Telescope Science Institute (STScI) which reproduces key functionality of the official GALEX pipeline in Python and makes it possible for individual researchers to generate the photon level data and calibrated lightcurves or integrated images. It also opens the possibility of modifying or improving upon the astrometric and photometric calibrations.

Additionally, the authors and MAST have undertaken to reduce all available GALEX data the time-tagged and sky-projected photon level events in a publicly accessible database. The database currently contains all photons collected in the GALEX direct imaging mode up through General Release (GR) 6/7 of the mission, comprising over 1 trillion events with a total volume over 100Tb. In addition to the standalone calibration pipeline, _gPipeline_, the authors have created tools ("the database tools") for querying and working with output from the photon database. These include
* _gFind_, for searching the database for specific coverage.
* _gAperture_, for extracting photon lists and generating calibrated lightcurves of specific targets.
* _gMap_, for creating calibrated images and movies.

##Getting Started
If you are a new user of the GALEX data or have not familiarized yourself with the quirks and intricacies of the GALEX detectors and calibration, please read the official [Technical Documentation](http://www.galex.caltech.edu/researcher/techdocs.html) with particular attention paid to [Chapter 3 - Pipeline Overview - Imaging](http://www.galex.caltech.edu/researcher/techdoc-ch3.html). This will answer many common questions of new users such as
* What is the difference between an observation, a visit, and an eclipse?
* What is a "dither correction?"
* How is a _relative response map_ different from a _flat field_?
* What is a "stim?"

The primary mission calibration paper (henceforth "the calibration paper") is also enlightening and thorough:

Morrissey, Patrick, et al. "The calibration and data products of GALEX." The Astrophysical Journal Supplement Series 173.2 (2007): 682.

###Installation Instructions
The standalone tools are written exclusively in Python, a flexible and powerful interpreted programming and data exploration language that has been increasingly adopted in many fields of research and especially astronomy. The "installation" in this case refers to installing the required version of Python and the non-standard Python modules ("dependencies") called by the gPhoton software. For naive Python users, we suggest simply [downloading the Anaconda distribution](https://store.continuum.io/cshop/anaconda/) of Python which contains all of the required dependencies. For advanced users and developers, we suggest that you manage the dependencies yourself; the complete dependency list and suggested installation instructions are given under **Manual Package Management** below.

Because the standalone tools are written in Python, they are theoretically cross platform. The software has only been thoroughly used and tested on Ubuntu Linux, however. At last attempt, we were not able to run the tools on Debian Linux because some of the required libraries were not yet supported.

####Anaconda + PIP (recommended for most users)
There are several versions of Python available which include not only the core "standard" version itself, but many common and popular modules as a single package, eliminating the need for users to manage such dependencies themselves. At present, the most promising of these appears to be Anaconda, which is available as a free download (with some advanced features available as paid add-ons). Anaconda contains all of the required dependencies for the gPhoton project. Advanced users and Python developers will probably want to manage their own dependencies; if you are not such a user, you can [download Anaconda here](https://store.continuum.io/cshop/anaconda/).

You must then install the most recent version of gPhoton from the Python Package Index (PyPI). Use the following terminal command: `pip install -i https://pypi.python.org/pypi gPhoton`

You will be able to run `gFind`, `gAperture`, `gMap`, and `gPipeline` as scripts directly from your terminal as well as `import gPhoton` from within a Python interactive session. The examples in this guide assume that gPhoton was installed using _pip_.

####Manual Package Management (advanced users)
Advanced users and developers will want more control over their installation and direct access to the source code and git repository. Importantly, this method will not give you access to the command line scripts unless you add them manually to your /usr/bin (or equivalent) directory.

#####Obtaining the Source Code
**For developers:** Obtain the source code by cloning the master branch of the [gPhoton repository on Github](https://github.com/cmillion/gPhoton). Instructions for getting started using Github can be found [here](https://help.github.com/categories/54/articles), and instructions specifically for cloning repositories can be found [here](https://help.github.com/articles/which-remote-url-should-i-use#cloning-with-ssh). Once you've cloned the repository, it will be straightforward for your to update your local version when we make updates to the master version.

Note: If you want to run the command line scripts (`gMap`, `gAperture`, `gFind`, and `gPipeline`) from your github checkout, you'll need to move them into the correct directory relative to the main module. They currently reside in the _bin_ directory directly under the repo (from which the PyPI installation puts them into _/usr/bin/_). You should copy them to the top level directory (which contains _gPhoton_, _docs_, etc.) and run them as scripts from there. (Running the scripts straight from _bin_ will result in an error like "ValueError: Attempted relative import beyond toplevel package")

#####Managing dependencies
You will need to install _python2.7_, _numpy_, _scipy_, _astropy_, _requests_ and _pandas_. The recommended commands for doing this appear below under the appropriate operating system.

The best specific tools for package installation and management shift rapidly. We'll try to keep this section up to date. If anything suggested here is actually _broken_, please let us know.

######Linux
Here are the recommended commands for Ubuntu. If you are using Fedora, substitute `yum` for `apt-get` everywhere.

    sudo apt-get install python-setuptools
    sudo apt-get install python-numpy python-scipy

You should use `pip` to get the latest versions of _requests_ and _astropy_. If _requests_ or _astropy_ is already installed, upgrade it by appending the `--upgrade` flag to the following calls.

    sudo pip install requests
    sudo pip install astropy
    sudo pip install pandas

######Mac (OSX)
**Draft.** For installing and managing your custom python build in Mac OSX, we suggest using the [MacPorts package](https://www.macports.org/). There is also a tutorial for installing Python on Mac with MacPorts [here](https://astrofrog.github.io/macports-python/).

    sudo port install py27-numpy
    sudo port install py27-scipy
    sudo port install py27-astropy
    sudo port install py27-requests
    sudo port install py27-pandas

Note: If your installation of requests complains about missing library dependencies, you may need to install them explicitly with the following command: `sudo pip install --upgrade pyopenssl ndg-httpsclient pyasn1`

######Windows
**Draft.** We haven't actually tried to do any of this on Windows. We suggest trying the [Enthought Python Distribution (EPD)](https://www.enthought.com/products/epd/) or the aforementioned [Anaconda](https://store.continuum.io/cshop/anaconda/) distribution.

##The Database Tools
The "database tools" or "photon tools" are the command line programs that provide basic functionality for interacting with he photon database at MAST to produce scientifically useful data products.

_Note:_ For the rest of this User Guide, we will use the M dwarf flare star GJ 3685A as our standard example target. The GALEX observation of this flare was described in _Robinson, et al. "GALEX observations of an energetic ultraviolet flare on the dM4e star GJ 3685A." The Astrophysical Journal 633.1 (2005): 447._ It's a good test because it has an obvious and dramatic light curve; you'll know it when you see it.

###gFind.py
_gFind_ is the data location tool. Given a target sky position (and, optionally, bands and time ranges), it will return the estimated raw exposure time and approximate time ranges of data that are currently available in the photon database. That is, _gFind_ is your convenient utility for assessing what data is currently available for use by _gAperture_ and _gMap_. Attempt the following command.

    gFind -b 'NUV' -r 176.919525856024 -d 0.255696872807351

This should report 2802.0 seconds of "raw" exposure time in five exposures along with the time ranges (in "GALEX Time") of those observations.

####Tweaking Search Parameters
The allowable gap between two time stamps for them to be considered "contiguous" is adjustable with the `--maxgap` parameter. Because MIS-depth GALEX observations are ~1600 seconds long, a good way to get time ranges that approximate the mission bookkeping of "visits" would be to set `--maxgap 1600`, and this is the default for this parameter in _gFind_. If you wanted only time ranges with no gaps, then a reasonable setting would be for one second,  because that is the default spacing between adjacent entries in the aspect table. As an example, try the following command to consider data with gaps of as much as 100 seconds to be contiguous.

    gFind -b 'NUV' -r 176.919525856024 -d 0.255696872807351 --maxgap 100

If, additionally, you want to exclude contiguous time ranges below some minimum raw exposure time, pass that time to the `--minexp` parameter.

    gFind -b 'NUV' -r 176.919525856024 -d 0.255696872807351 --minexp 200

And, naturally, the `--gap` and `--minexp` parameters can be used in conjunction.

    gFind -b 'NUV' -r 176.919525856024 -d 0.255696872807351 --maxgap 100 --minexp 200

If you want to exclude times when the source is on the edge of the detector, you can adjust the `--detsize` parameter to the desired effective detector _width_ (default = 1.1 degrees, to exclude the edges with a wide buffer).

    gFind -b 'NUV' -r 176.919525856024 -d 0.255696872807351 --detsize 0.5

Note that this reduces the depth of some exposures in addition to eliminating some exposures entirely. This makes sense in light of the fact that the position of a GALEX source on the detector moves during an observation.

_For the curious:_ The estimates returned by _gFind_ are computed by finding the boresight pointings (in the aspect database) which fall within a detector _radius_ (nominally 0.625 degrees) of the desired sky position and comparing the associated time stamps against the time stamps of data that has actually been loaded into the photon database. The _predicted_ keyword performs this same search on the aspect solutions only without comparing it against the database.

####Alternative I/O Formats
Rather than passing RA (`-r`) and Dec (`-d`) separately, you can pass them to `--skypos` as follows.

    gFind -b 'NUV' --skypos '[176.919525856024,0.255696872807351]'

If you are interested only in the available raw exposure times and not the individual time ranges, pass the `--exponly` flag.

    gFind -b 'NUV' --skypos '[176.919525856024,0.255696872807351]' --exponly

####Calling from within the Python Interpreter
If you're creating scripts using the gPhoton tools, it might be more convenient to work from directly within the Python interpreter. _gFind_ can be imported like any other module in the following manner.

    import gPhoton.gFind

For information about the module:

    help(gPhoton.gFind)

Note that the heavy lifting within gFind is performed by a function _itself_ called gFind. For example, an equivalent call to

    gFind -b 'NUV' -r 176.919525856024 -d 0.255696872807351 --gap 100 --minexp 100

within the interpreter would be the following.

    import gPhoton.gFind
    gPhoton.gFind.gFind(band='NUV',skypos=[176.919525856024,0.255696872807351],maxgap=100.,minexp=100.)

###gAperture.py
_gAperture_ is the photometry tool which computes source fluxes or light curves for specified targets and time ranges with customizable apertures and background annuli. If an output filename is provided, the light curve data will be written to a .csv file.

The minimum required parameters are RA (`-r` or `--ra`), Dec (`-d` or `--dec`), and aperture radius (`-a`), all in decimal degrees. This will compute the integrated flux over all available data with no background subtraction. For our flare star example and an aperture with radius of 0.03 degrees, that command looks like this.

    gAperture -b 'NUV' -r 176.919525856024 -d 0.255696872807351 -a 0.03

This will print out the integrated magnitudes for the five observations listed by gFind above. Note a ~3 magnitude jump in brightness in the third (middle) observation; that's the flare.

You can limit the computation to specific time ranges with `--t0` and `--t1` or `--trange` (or `--tranges`). You can also perform a background correction by specifying the inner and outer radii of a background annulus (in decimal degrees) centered on the target with the `-i`/`--inner` and `-o`/`--outer` or `--annulus` flags.

Let's limit further analysis to the time range in which flare occurs in this data, `[766525332.995,766526576.995]`. Let's also extract the background from an annulus with an inner radius of 0.03 degrees and an outer radius of 0.04 degrees.

    gAperture -b 'NUV' -r 176.919525856024 -d 0.255696872807351 -a 0.03 -i 0.03 -o 0.04 --t0 766525332.995 --t1 766526576.995

which is equivalent to

    gAperture -b 'NUV' --skypos [176.919525856024,0.255696872807351] -a 0.03 --annulus [0.03,0.04] --trange [766525332.995,766526576.995]

Rather than simply printing magnitudes to the terminal, you can pass a .csv filename to gAperture with the `--filename` or `-f` flag, into which a significant amount of data about the observation will be written. (See column definitions below.)

    gAperture -b 'NUV' -r 176.919525856024 -d 0.255696872807351 -a 0.03 -i 0.03 -o 0.04 --t0 766525332.995 --t1 766526576.995 -f 'lightcurve.csv'

Note that if you try to run that command a second time, it won't let you because the software detects that the file already exists and suggests that you use `--clobber` (or `--overwrite`) to force it to overwrite _lightcurve.csv_. Therefore, `--clobber` is appended to all of the following commands to avoid this error. (But you should use it with caution. It's there for a reason.)

If you want to generate a light curve rather than an integrated value, pass the desired (raw) bin depth in seconds to the step size flag (`-s`). For example, to generate a light curve with 100 second bins, try the following.

    gAperture -b 'NUV' -r 176.919525856024 -d 0.255696872807351 -a 0.03 -i 0.03 -o 0.04 --t0 766525332.995 --t1 766526576.995 -f 'lightcurve.csv' -s 100. --overwrite

For any command, you can always request more information be printed to the terminal by setting the `--verbose` or `-v` flag to a number between 1-3 (defualt is 0) where larger numbers indicate increasing levels of output. Setting `-v 3` will print out complete SQL commands and should really only be used for debugging.

####Lightcurve File Column Definitions
**NOTE:** The column definitions for the .csv output from _gAperture_ are in flux. These are the column definitions as of the v1.27.0 build.

**NOTE:** The columns are not necessarily written to the output file in the order given. Nor are the columns necessarily fixed in order at all. You should parse the lightcurve file on the column _name_ and not the column number.

**NOTE:** Each row of the light curve file corresponds to a single time bin. So assume that every description below includes the implicit appendix of "within the time bin" unless stated otherwise.

**NOTE:** All times ware in GALEX seconds. Positions are in degrees. Areas are in square degrees. Fluxes are in units of erg sec^-1 cm^-2 Å^-1.

1. flat_counts - The sum over all flat-corrected counts within the aperture.
2. mcat_bg - Estimated background brightness as pulled from the visit-level MCAT and scaled to the area of the aperture.
3. bg_counts - Raw number of counts within the background annulus.
4. flux_bgsub_err - Estimated 1-sigma error in `flux_bgsub` value.
5. cps_mcatbgsub - Countrate within the aperture, corrected for background using the visit-level MCAT background estimates.
6. counts - Total number of uncorrected counts within the photometric aperture.
7. mag_mcatbgsub - AB Magnitude within the aperture, corrected for background using the visit-level MCAT background estimates.
8. cps_err - Estimated background in the countrate within the aperture, assuming no contribution from background. (i.e. sqrt(n))
9. mag_bgsub - AB Magnitude within the aperture, corrected by the background estimated from the annulus.
10. cps_bgsub - Countrate within the aperture, corrected by the background estimated from the annulus.
11. detys - Mean detector Y position of events within the aperture.
12. flux_bgsub - Flux within the aperture, corrected by the background estimated from the annulus.
13. flux_err - Estimated 1-sigma error in `flux` value.
14. mag_err_1 - Estimated upper 1-sigma error in `mag` value.
15. cps_bgsub_err - Estimated 1-sigma error in `cps_bgsub` value.
16. t1_data - Final timestamp of events within the aperture.
17. bg - Contribution of background, as estimated from the annulus and scaled to the area of the aperture.
18. responses - Mean value of the flat assigned to events within the aperture.
19. t_mean - Mean timestamp of events within the aperture.
20. cps_mcatbgsub_err - Estimated 1-sigma error on `cps_mcatbgsub` value.
21. mag_bgsub_err_1 - Estimated upper 1-sigma error on `mag_bgsub` value.
22. mag_err_2 - Estimated lower 1-sigma error on `mag` value.
23. t0_data - Earliest timestamp of events within the aperture.
24. racent - Mean right ascension of events within the aperture.
25. deccent - Mean declination of events within the aperture.
26. mag - AB Magnitude within the aperture, uncorrected for background.
27. exptime - Estimated effective exposure time (correct for dead time and shutter).
28. bg_flat_counts - Total of flat-corrected counts within the background annulus.
29. detxs - Mean detector X position of all events within the aperture.
30. t0 - Lower time delimiting the bin.
31. t1 - Upper time delimiting the bin.
32. mag_mcatbgsub_err_2 - Estimated lower 1-sigma error on `mag_mcatbgsub` value.
33. flux - Flux witin the aperture, uncorrected for background.
34. mag_mcatbgsub_err_1 - Estimated upper 1-sigma error on `mag_mcatbgsub` value.
35. flags - Automatically generated gAperture quality flag. Bins with a flag that is non-zero should not be naively trusted. See flag definitions below for more information.
36. mag_bgsub_err_2 - Estimated lower 1-sigma error on `mag_bgsub`.
37. detrad - Mean detector radius (distance from detector center) for events within the aperture.
38. cps - Countrate within the aperture, uncorrected for background.
39. flux_mcatbgsub_err - Estimated 1-sigma error on `flux_mcatbgsub` value.
40. flux_mcatbgsub - Flux within the aperture, corrected for background using the visit-level MCAT values.

#####Flag Column Definitions
These flags are automatically set in software based upon conditions that we know to reproducibly generate misleading lightcurves. The flags are additive in binary, so it's possible to have more than one flag set at a time. They are defined as follows:

1 - 'hotspot' - events in pixels contiguous to a hotspot masked region

2 - 'mask edge' - events in pixels contiguous to the detector edge

4 - 'exptime' - bin contains < 50% exposure time coverage

8 - 'respose' - events weighted with response < 0.7

16 - 'nonlinearity' - local countrate exceeds 10% response dropoff

32 - 'detector edge' - events outside of 0.5 degrees of detector center

####Calling from within the Python Interpreter
You can also import and work with _gAperture_ and its modules from within the Python interpeter.

    import gPhoton.gAperture

###gMap.py

**NOTE: As of v1.26.0, gMap does not apply any exposure time correction. The images (including intensity) should be used for assessment purposes, but not photometric analysis.**

_gMap_ is the image creation tool. It can generate integrated count, intensity, and response (equivalent to GALEX _cnt_, _int_ and _rrhr_) maps of arbitrary size<sup>+</sup>, shape and depth, including coadds across epochs and survey designation. It can also create "movie" (time-binned, multi-plane) versions of such maps. _(Note that in v1.23.0, image creation is quite slow because it uses an older method of interpolating the response map from the flat rather than the new and faster method used by gAperture of weighting the photons directly by the flat. This will be fixed in v1.24.0)_

_gMap_ writes all image files in the Flexible Image Transport System (FITS) standard, which is an uncompressed archival data format favored by many astronomy applications (and was used by the GALEX mission for archival products). FITS images generated by _gMap_ have headers which describe the World Coordinate System (WCS) information defining their orientation on the sky as well as effective exposure time and other metadata for the observation as a whole. For multi-frame images (movies), the per-plane information is described in a table in the FITS secondary HDU; this table describes start time, stop time, and effective exposure for each frame (within appropriately labelled columns).

<sup>+</sup>Caveat emptor. "Arbitrariness" is limited by your available patience and RAM.

####Count Maps
Count (_cnt_) maps are integrated, aspect corrected but uncalibrated (that is, not adjusted for resonse or exposure time) images of the sky. Count images are good for "quick looks" at the data, to ensure that you are pointing in the location that you expected and that you are seeing the sources or features desired. But because they are not calibrated by either the relative response or the effective exposure time, you should not use them for photometric analyses of any kind.

You can create a count image from the command line by specifying the band (`-b`), sky position (`-r` and `-d` or just `--skypos`), the angular extent of the desired image in degrees (`--angle`), an output FITS filename (`--count`), and optionally a time range (`--t0` and `--t1` or just `--trange` or `--tranges`). Try the following simple command.

    gMap -b 'FUV' -r 176.919525856024 -d 0.255696872807351 --angle 0.5 --t0 766525332.995 --t1 766526576.995 --count 'count.fits'

Or, using the alternative formats for specifying sky position and time range, try the following.

    gMap -b 'FUV' --skypos '[176.919525856024,0.255696872807351]' --angle 0.5 --tranges '[[766525332.995, 766526576.995], [919755500.995, 919755916.995]]' --count 'count.fits'

Note that the most recent command created an image with two planes corresonding to the time ranges `[766525332.995, 766526576.995]` and `[919755500.995, 919755916.995]`, respectively. This behavior can be used to specify custom time ranges for a movie (see **Movies** below), or you can force _gMap_ to coadd these two time ranges with the `--coadd` flag.

    gMap -b 'FUV' --skypos '[176.919525856024,0.255696872807351]' --angle 0.5 --tranges '[[766525332.995, 766526576.995], [919755500.995, 919755916.995]]' --count 'count.fits' --coadd

However, if you do not specify any time range, _gMap_ will automatically use all available exposure. (And uses the same underlying function as _gFind_ to locate said exposure time, so therefore will also accept the `--maxgap` and `--minexp` and `--detsize` keywords if desired.)

    gMap -b 'FUV' --skypos '[176.919525856024,0.255696872807351]' --angle 0.5 --count 'count.fits' --coadd

####Intensity Maps
Intensity maps (_int_) are integrated, aspect corrected images which have been corrected for both relative response as well as effective exposure time. If you need to perform photometric analysis on an image, you need an intensity map. Note the warning above about long run times for generating the _response_ for this. If you want to create an intensity map, pass a FITS filename to the `--intensity` flag.

    gMap -b 'FUV' --skypos '[176.919525856024,0.255696872807351]' --angle 0.5 --intensity 'intensity.fits' --coadd

####Movies
You can turn any of the above map types into a _movie_ (that is, a time-binned, multi-plane image), by passing the desired (raw) bin depth in seconds to the `--frame` parameter. For example, to create a count image with 100 second depth image planes, try the following.

    gMap -b 'FUV' -r 176.919525856024 -d 0.255696872807351 --angle 0.5 --t0 766525332.995 --t1 766526576.995 --count 'count.fits' --frame 100

##gPipeline.py - The Standalone Calibration Pipeline
You will need to download the sample eclipse directory from [here](https://archive.stsci.edu/prepds/gphoton/cal/e31000.tar.gz). This directory contains the raw science (raw6), spacecraft state (scst), and refined aspect (asprta) files for eclipse e31000. Unzip this test eclipse into a known location, and then cd into that directory.

Run the first two example gPhoton commands below. These will generate all of the aspect corrected FUV and NUV photons for eclipse 31000 as comma separate value (.csv) files. While gPhoton is running, the terminal will update with the photon chunk and current processing rate, mostly just so you know that something is happening. The photons will be dumped into files named [NF]UVPhotons.csv which will be quite large (several Gb), so make sure that you have enough disk space available. It can also take minutes to hours (depending on the observation and speed of your computer); if processing is taking longer than you want, you can kill gPipeline at any time with `Ctrl+C`. New runs will overwrite .csv output from previous runs.

This generates the aspect corrected FUV photons for eclipse 31000 and writes them to `FUVphotons.csv`.

    gPipeline -r 'MISWZN11_12494_0315_0002-fd-raw6.fits' -a 'MISWZN11_12494_0315_0002-asprta.fits' -o 'FUVphotons' -s 'MISWZN11_12494_0315_0002-scst.fits' -d 'SSD_fuv_31000.tbl'

This generates the aspect corrected NUV photons for eclipse 31000 and writes them to `NUVphotons.csv`. Note that because GALEX countrates are much higher in the NUV, this will take ~10x longer to run than the previous command.

    gPipeline -r 'MISWZN11_12494_0315_0002-nd-raw6.fits' -a 'MISWZN11_12494_0315_0002-asprta.fits' -o 'NUVphotons' -s 'MISWZN11_12494_0315_0002-scst.fits' -d 'SSD_nuv_31000.tbl'

###Optional parameters
A number of the command line parameters given above are actually option. If (and only if) you have a working internet connection, then the aspect file (`-a`) parameter can be omitted; the software will instead query the aspect database table at MAST for the appropriate values. The Stim Separation Data (SSD) file parameter (`-d`) is _always optional_ because the values in this reference table can be generated directly from the raw data (at a very small cost in terms of total run time). Try the following command, which omits both of these parameters.

    gPipeline -r 'MISWZN11_12494_0315_0002-fd-raw6.fits' -o 'FUVphotons' -s 'MISWZN11_12494_0315_0002-scst.fits'

###NULL vs. non-NULL data
A large number of events cannot be aspect corrected either because they fall in a time range for which noa spect solution is available (for a variety of possible reasons) or because they don't actually fall on the detector proper (e.g. stims). Such events appear in the .csv output file with _NULL_ entries for RA and Dec and are referred to as "null data." These data are of limited scientific interest and it was convenient for us to separate them in the MAST database, so you can optionally write null data to a separate file from non-null data by passign the `-u` flag. This will create a second .csv file for null data with "_NULL.csv" appended to the filname.

    gPipeline -r 'MISWZN11_12494_0315_0002-fd-raw6.fits' -a 'MISWZN11_12494_0315_0002-asprta.fits' -o 'FUVphotons' -s 'MISWZN11_12494_0315_0002-scst.fits' -u

###Multi-visit Eclipses
For multi-visit eclipses (e.g. AIS), you can specify more than one aspect file for a single raw6 file using the following syntax.

    -a '../FOO_sv01-asprta.fits,../FOO_sv02-asprta.fits,../FOO_sv03-asprta.fits'

Again, if you do not specify `-a` at all, then the software will query the aspect database at MAST for the appropriate values, even for multi-visit eclipses.

###Other Notes
####Photon File Column Definitions
The column definitions for the .csv file output by gPhoton.py are as follows.

1. _t_ - time of the event (in "GALEX Time" = "UNIX Time" - 315964800)
2. _x_ - detector x position
3. _y_ - detector y position
4. _xa_ - (if you don't already know, it's not important)
5. _ya_ - (ditto)
6. _q_ - (also ditto)
7. _xi_ - aspect corrected detector position
8. _eta_ - aspect corrected detector position
9. _ra_ - aspect corrected right ascension of the event (decimal degrees)
10. _dec_ - aspect corrected declination of the event (decimal degrees)
11. _flags_ - status information (see below)

####Flag Column Definitions
These are the definitions of various values of the _flag_ column in the gPhoton .csv output file. Note that these are not identical to any flag column definitions for output created by the canonical pipeline. In general, end users will be most interested in events for which _flag = 1_, and this is the default search criterion for all queries performed by the Photon Tools.

0. successfully calibrated (no errors)
1. _N/A_
2. event skipped
3. _N/A_
4. _N/A_
5. classified as stim
6. hotspot (covered by hotspot mask)
7. bad aspect (unkown aspect solution)
8. out of range (off of detector prior to the wiggle correction)
9. bad walk (off of detector prior to the walk correction)
10. bad linearity (off of the detector prior to the linearity correction)
11. bad distortion (off of the detector prior to the distortion correction)
12. aspect jump (questionable aspect due to occuring during a jump or gap in the aspect solution or bracketed by one or more flagged aspect values)
13. _N/A_
14. _N/A_
15. _N/A_

##Common Questions, Issues, and Gotchas
1. **"My data is not available!"** You can verify that data for your desired target does or does not exist in the database and present by using the `gFind` commands described above. If data for your target is not available, there are two possible explanations: (1) we have not yet loaded those observations into the photon database, or (2) that target was never observed by the GALEX mission. As of this writing, we have loaded all direct imaging data up through the GR6/7 release, which corresponds to the end of the NASA-funded mission and beginning of the CAUSE phase. You can confirm that your target was, indeed, observed by GALEX by searching for it in the [GALEX Catalog](http://galex.stsci.edu/GR6/?page=mastform).
2. **Exposure Time.** (Or **"Why isn't the exposure time equal to the bin width?"**) When calculating the exposure time of a GALEX observation, one cannot simply subtract the end time from the start time. This _raw exposure time_ (t<sub>raw</sub>) must be adjusted into an _effective exposure time_ (t<sub>eff</sub>) which takes into account certain microchannel plate detector properties. In particular, the raw exposure time must be adjusted by both the _shutter_ (t<sub>shut</sub>) and the _deadtime ratio_ such that: t<sub>eff</sub> = (t<sub>raw</sub> - t<sub>shut</sub>) x deadtime. In general, you can expect t<sub>eff</sub> ~= 0.8 * t<sub>raw</sub>, though the actual ratio varies wildly as a function of field brightness and the intersection of bin boundary times with observation time ranges.
    * **Shutter Correction.** GALEX did not observe all parts of the sky at all times. Even when GALEX was observing a particular part of the sky, there were times when little or no usable data was actually recorded. Data for specific targets or time ranges might exist but still be (practically) unusable for a variety for a variety of reasons which may include (1) that the spacecraft was not observing at nominal high voltage ("HVNOM"), (2) the spacecraft was observing a different part of the sky or in spectroscopic ("grism") mode, (3) a valid aspect solution could not be reconstructed from the available data (perhaps due to too few detectable reference stars), or (4) there was a temporary gap in the data due to a spacecraft or data transmission
    anomaly. The _shutter correction_ accounts for these gaps by conceptualizing a "virtual shutter" that is considered closed whenever a gap of 0.05 seconds occurs in the data covering a particular source (as determined by boresight center and detector width) and within a particular time range. The integrated time of all such 0.05 second gaps in any given time range (t<sub>shut</sub>) is subtracted from t<sub>raw</sub>.
    * **Deadtime Correction** The GALEX microchannel plates could only detect a single event at a time; that is, when one photon event was being read into the electronics, any othe rincident photon events were missed entirely. This effectively reduced the exposure time of any observation by a ratio known as the "deadtime correction" which scaled as a function of field brightness (i.e. global count rate). Formally, the deadtime correction is the estimated fraction of the raw exposure time that the detector is not properly observing because it is engaged in readout. The deadtime changes from observation to observation and even moment to moment and, so, is ideally calculated continuously. There are two ways to estimate the deadtime.
        1. **Direct Measurement.**
        2. **Empirical Estimate.**
3. **Relative Response.** Microchannel plates do not have true "pixels" in the same sense as CCD or CMOS detectors. Microchannel plates are also prone to "burn in" or _local gain sag_ where regions of the detector exposured to high countrate sources can be local depleted of charge and suffer either temporary or permanent loss of sensitivity. To avoid this local gain sag, the telescope boresight did not stare fixedly at a point on the sky over the course of any observation, and this relatively motion was then accounted for in the _aspect correction_ stage of the calibration pipeline. As such, to apply the detector flat field in sky-coordinates, it must be interpolated along the effective field of view during the observation to form a _relative response_ map. **[...]**
4. **Gain Sag.** [COMING SOON]
5. **What are the photometric errors of the data?** To good approximation, the photometric errors are equal to sqrt(counts)/t<sub>eff</sub>, which is just the Gaussian error on the point source. Technically, one should include the background error in quadrature, but if the sky background is large enough that it's contributing significantly to your error, then something has gone wrong with the background correction or your source has a low enough signal-to-noise that you should be performing some more careful error analysis. Note that developing a more thorough end-to-end error model is on our todo list, but near the bottom.
6. **How is background estimation performed?** The counts per second per unit area are computed for the background annulus after masking out known catalog sources (down to mag 22 as the default) and then scaling that to the area of the source aperture. In the future, we may support use of the catalog background values, but this won't correct for the residual zodiacal light which can change over the course of an observation and potentially produce a false signal of variability.
7. **Why are there negative counts per second or NaN magnitudes?** This happens when the measured background (in counts per second per area) is larger than the measured source (in counts per second per area) such that `cps<sub>source</sub>-cps<sub>background</sub><0`. In most cases, this is for the obvious reason that the signal is near or has dropped below the background. It could also be due to a failure of the background estimation algorithm to mask out sources in the annulus.
8. **How do I convert GALEX time stamps to something meaningful?** The GALEX time stamps are in "GALEX Time" which is is defined UNIX Time less 315964800 seconds. (t<sub>GALEX</sub> = t<sub>UNIX</sub> - 315964800) UNIX Time is a standard defined as the number of seconds that have elapsed since January 1, 1970. A number of utilities exist online for converting UNIX time to something meaningful (like Julian Date).
