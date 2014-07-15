#DRAFT DRAFT DRAFT DRAFT DRAFT

**User Documentation for the MAST/GALEX Photon Database and Tools**

Chase Million<sup>1</sup>, Bernie Shiao<sup>2</sup>, Scott Fleming<sup>2</sup>,
Myron Smith<sup>2</sup>

<sup>1</sup>Million Concepts (chase.million@gmail.com),
<sup>2</sup>Space Telescope Science Institute

###Raison d'Etre

The MAST/GALEX photon database and tools exist in an effort to maximize the flexibility and utility of the GALEX data set. The GALEX detectors are microchannel plates which record detector position and time-of-arrival information for every detected photon event with a time resolution of five thousandths of a second, composing a huge and rich short time domain survey of the UV. Due to digital storage space and processing limitations, the data was only formally released by the mission team as integrated images. The photon list files--internally known as _x-files_--were only available by special request, and there was little to no additional support for their calibration or use.

The official GALEX calibration pipeline software ("the canonical pipeline")--written in about half a dozen langauges with a sprawling network of complex dependencies--has also never been successfully ported to any system outside of the GALEX internal network at Caltech. Even though the source code for this pipeline will be made publicly available through the Mikulski Archive at Space Telescope (MAST) in the future, the end of the GALEX project would have effectively marked the end of the capability to generate photon level data specifically and revisit the GALEX calibration more generally. A software tool known as _gPhoton_ ("the standalone pipelien") has been developed by the authors with support of MAST and Space Telescope Science Institute (STScI) which reproduces a large portion of the official GALEX pipeline in Python and makes it possible for individual researchers to both generate their own photon level data and calibrated lightcurves or integrated images, and also to adjust or improve upon the calibration.

Additionally, the authors and MAST have undertaken to process all of the GALEX data with gPhoton and store the photon level data in a database. Once the database is fully populated (_est. 2015_), it will comprise ~180 terabytes and contain approximately 1.5 trillion rows (at one event per row). In addition to the standalone calibration pipeline, _gPhoton_, the authors have created tools ("the database tools") for querying and working with output from the photon database. These include
* _gFind_, for searching the database for specific coverage.
* _gAperture_, for extracting photon lists and generating calibrated lightcurves of specific targets.
* _gMap_, for creating calibrated images and movies.

##Getting Started
If you are a new user of the GALEX data or have not familiarized yourself with the quirks and intricacies of the GALEX detectors and calibration, please read the official [Technical Documentation](http://www.galex.caltech.edu/researcher/techdocs.html). Please pay particular attention to [Chapter 3 - Pipeline Overview - Imaging](http://www.galex.caltech.edu/researcher/techdoc-ch3.html). This will answer many common questions of new users such as
* What is the difference between an observation, a visit, and an eclipse?
* What is a "dither correction?"
* How is a _relative response map_ different from a _flat field_?
* What is a "stim?"

###High Level Description of the Software
**[PLACEHOLDER TEXT]**

###Obtaining the Source Code
**Preferred method:** Obtain the source code by cloning the master branch of the [gPhoton repository on Github](https://github.com/cmillion/gPhoton). Instructions for getting started using Github can be found [here](https://help.github.com/categories/54/articles), and instructions specifically for cloning repositories can be found [here](https://help.github.com/articles/which-remote-url-should-i-use#cloning-with-ssh). Once you've cloned the repository, it will be straightforward for your to update your local version when we make updates to the master version.

**Alternative method:** If you can't set up Github on your machine for some reason (such as permissions issues), you can download the master branch with the [**Download Zip**](https://github.com/cmillion/gPhoton/archive/master.zip) button on Github. This will initiate a large (~200Mb) download of a file called `master.zip` which will contain the source code. Note that in order to obtain software updates, you will need to re-download and re-extract this file.

###Installation Instructions
The standalone tools are written exclusively in Python which is a flexible and powerful interpreted programming and data exploration language that is being increasingly adopted in many fields of research and especially astronomy. The "installation" in this case refers to installing the required version of Python and the non-standard Python modules ("dependencies") called by the gPhoton software. For naive Python users, we suggest simply [downloading the Anaconda distribution](https://store.continuum.io/cshop/anaconda/) of Python which contains all of the required dependencies built-in; Anaconda is described in greater detail under **Prebuilt Distributions** below. For advanced users and developers, we suggest that you manage the dependencies yourself; the complete dependency list and suggested installation instructions are given under **Manual Package Management** below.

Because the standalone tools are written in Python, they are theoretically cross platform. The software has only been thoroughly used and tested on Ubuntu Linux, however. At last attempt, we were not able to run the tools on Debian Linux because some of the required libraries are not yet supported.

####Prebuilt Distributions
There are several versions of Python available which include not only the core "standard" version itself, but many common and popular modules as a single package, eliminating the ened for users to manage such dependencies themeselves. At present, the most promising of these appears to be Anaconda, which is available as a free download (with some advanced features available as paid add-ons). Anaconda contains all of the required dependencies for the gPhoton project. Advanced users and Python developers will probably want to manage their own dependencies; if you are not one of these, you can [download Anaconda here](https://store.continuum.io/cshop/anaconda/) and get started using gPhoton immediately.

####Manual Package Management
You will need to install _python2.7_, _numpy_, _scipy_, _astropy_ and _requests_. To use some of the Python helper utilities (contained in gphoton_utils.py), you will also need to install _pandas_. The recommended commands for doing this appear below under the appropriate operating system.

The best specific tools for package installation and management shift rapidly. We'll try to keep this section up to date. If anything suggested here is actually _broken_, please let us know.

#####Linux
Here are the recommended commands for Ubuntu. If you are using Fedora, substitute `yum` for `apt-get` everywhere.

    sudo apt-get install python-setuptools
    sudo apt-get install python-numpy python-scipy

You should use `pip` to get the latest versions of _requests_ and _astropy_. If _requests_ or _astropy_ is already installed, upgrade it by appending the `--upgrade` flag to the following calls.

    sudo pip install requests
    sudo pip install astropy
    sudo pip install pandas

#####Mac (OSX)
**DRAFT** For installing and managing your custom python build in Mac OSX, we suggest using the [MacPorts package](https://www.macports.org/). There is also a tutorial for installing Python on Mac with MacPorts [here](https://astrofrog.github.io/macports-python/).

    sudo port install py27-numpy
    sudo port install py27-scipy
    sudo port install py27-astropy
    sudo port install py27-requests
    sudo port install py27-pandas

#####Windows
**[PLACEHOLDER TEXT]** We haven't actually tried to do any of this on Windows.

###Testing Your Build
If you want to test your build or run any of the `gPhoton` commands below, you will need to download the sample eclipse directory from [here](https://www.dropbox.com/s/2c26jafccqz5ahh/e31000.tar.gz). This directory contains the raw science (raw6), spacecraft state (scst), and refined aspect (asprta) files for eclipse e31000. Unzip this test eclipse into the same directory as gPhoton (i.e. the directory `e31000` should be on the same level and in the same directory as `source` and `cal`).

From the command line, navigate to the `gPhoton/source` directory. Run the first two example gPhoton commands below. These will either generate errors (if a package is missing or wrong) or they will generate all of the aspect corrected FUV and NUV photons for eclipse 31000 as comma separate value (CSV for .csv) files. While gPhoton is running, the terminal will update with the photon chunk and current processing rate, mostly just so you know that something is happening. The photons will be dumped into files named [NF]UVPhotons.csv which will be quite large (several Gb), so make sure that you have enough disk space available. If a few chunks are processed without errors, your installation is likely fine. You can kill gPhoton at any time with `Ctrl+C` New runs will overwrite .csv output from previous runs.

**Convenient hack:** The .csv files get updated frequently, so you can generate a "sample" .csv file by letting gPhoton run for a few minutes and then killing it.

**Note:** If you don't have `PYTHONPATH` defined, then you will need to put `python` in front of all of the command line tool calls in this document. (e.g. `python ./gPhoton.py [...] [...] [...]`).

From the command line, navigate to the `gPhoton/source` directory. Then try running the following commands.

This generates the aspect corrected FUV photons for eclipse 31000 and writes them to `FUVphotons.csv`.

    ./gPhoton.py -r '../e31000/MISWZN11_12494_0315_0002-fd-raw6.fits' -a '../e31000/MISWZN11_12494_0315_0002-asprta.fits' -c '../cal/' -o '../e31000/FUVphotons' -s '../e31000/MISWZN11_12494_0315_0002-scst.fits' -d '../e31000/SSD_fuv_31000.tbl'

This generates the aspect corrected NUV photons for eclipse 31000 and writes them to `NUVphotons.csv`. Note that because GALEX countrates are much higher in the NUV, this will take ~10x longer to run than the previous command.

    ./gPhoton.py -r '../e31000/MISWZN11_12494_0315_0002-nd-raw6.fits' -a '../e31000/MISWZN11_12494_0315_0002-asprta.fits' -c '../cal/' -o '../e31000/NUVphotons' -s '../e31000/MISWZN11_12494_0315_0002-scst.fits' -d '../e31000/SSD_nuv_31000.tbl'

##gPhoton.py - The Standalone Calibration Pipeline

_Syntax Note:_ The name _gPhoton_ can refer to both the standalone calibration pipeline  `gPhoton.py` and the MAST/GALEX photon database project as a whole. We know that this is confusing. For this reason, we try to refer to `gPhoton.py` as "the standalone calibration pipeline" and reserve "gPhoton" to refer to the project as a whole.

###Optional parameters
A number of the command line parameters given above are actually option. If (and only if) you have a working internet connection, then the aspect file (`-a`) parameter can be omitted; the software will instead query the aspect database table at MAST for the appropriate values. The Stime Separation Data (SSD) file parameter (`-d`) is _always optional_ because the values in this reference table can be generated directly from the raw data (at a very small cost in terms of total run time). Try the following command, which omits both of these parameters.

    ./gPhoton.py -r '../e31000/MISWZN11_12494_0315_0002-fd-raw6.fits' -c '../cal/' -o '../e31000/FUVphotons' -s '../e31000/MISWZN11_12494_0315_0002-scst.fits'

###NULL vs. non-NULL data
A large number of events cannot be aspect corrected either because they fall in a time range for which noa spect solution is available (for a variety of possible reasons) or because they don't actually fall on the detector proper (e.g. stims). Such events appear in the .csv output file with _NULL_ entries for RA and Dec and are referred to as "null data." You can optionally write null data to a separate file from non-null data by passign the `-u` flag. This will create a second .csv file for null data with "_NULL.csv" appended to the filname.

    ./gPhoton.py -r '../e31000/MISWZN11_12494_0315_0002-fd-raw6.fits' -a '../e31000/MISWZN11_12494_0315_0002-asprta.fits' -c '../cal/' -o '../e31000/FUVphotons' -s '../e31000/MISWZN11_12494_0315_0002-scst.fits' -u

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

1. successfully calibrated (no errors)
2. _N/A_
3. event skipped
4. _N/A_
5. _N/A_
6. classified as stim
7. hotspot (covered by hotspot mask)
8. bad aspect (unkown aspect solution)
9. out of range (off of detector prior to the wiggle correction)
10. bad walk (off of detector prior to the walk correction)
11. bad linearity (off of the detector prior to the linearity correction)
12. bad distortion (off of the detector prior to the distortion correction)
13. aspect jump (questionable aspect due to occuring during a jump or gap in the aspect solution or bracketed by one or more flagged aspect values)
14. _N/A_
15. _N/A_
16. _N/A_

##The Database Tools
The "database tools" or "photon tools" are the command line programs that provide basic functionality for interacting with he photon database at MAST to produce scientifically useful data products.

_Note:_ For the rest of this User Guide, we will use the M dwarf flare star GJ 3685A as our standard example target. The GALEX observation of this flare was described in _Robinson, et al. "GALEX observations of an energetic ultraviolet flare on the dM4e star GJ 3685A." The Astrophysical Journal 633.1 (2005): 447._ It's a good test because it has an obvious and dramatic light curve; you'll know it when you see it.

###gFind.py
_gFind_ is the data location tool. Given a target sky position (and, optionally, bands and time ranges), it will return the estimated raw exposure time and approximate time ranges of data that is currently available in the photon database. That is, _gFind_ is your convenient utility for assessing what data is currently available for use by _gAperture_ and _gMap_. Attempt the following command.

    ./gFind.py -b 'NUV' -r 176.919525856024 -d 0.255696872807351

_For the curious:_ The estimates returned by _gFind_ are computed by finding the boresight pointings (in the aspect database) which fall within a detector radius (~0.625 degrees) of the desired sky position and comparing the associated time stamps against the time stamps of data that has actually been loaded into the photon database.

####Tweaking Search Parameters
The default allowable gap between two time stamps for them to be considered contiguous for the purposes of a _gFind_ estimate is one second. This is a reasonable minimum gap because it is the default spacing between adjacent entries in the aspect table. This parameter is adjustable, however, with the `-g` or `--gap` parameter. To consider data with gaps of as much as 100 seconds to be contiguous--a reasonable value if you want to treat all data in the same eclipse as part of the same observation--try the following command.

    ./gFind.py -b 'NUV' -r 176.919525856024 -d 0.255696872807351 --gap 100

If, additionally, you want to exclude contiguous time ranges below some minimum raw exposure time, pass that time to the `--minexp` parameter.

    ./gFind.py -b 'NUV' -r 176.919525856024 -d 0.255696872807351 --minexp 100

And, naturally, the `--gap` and `--minexp` parameters can be used in conjunction.

    ./gFind.py -b 'NUV' -r 176.919525856024 -d 0.255696872807351 --gap 100 --minexp 100

If you want to exclude times when the source is on the edge of the detector, you can adjust the `--detsize` parameter to the desired effective detector _width_ (default = 1.25 degrees).

    ./gFind.py -b 'NUV' -r 176.919525856024 -d 0.255696872807351 --detsize 1

_For the curious:_ The estimates returned by _gFind_ are computed by finding the boresight pointings (in the aspect database) which fall within a detector _radius_ (nominally 0.625 degrees) of the desired sky position and comparing the associated time stamps against the time stamps of data that has actually been loaded into the photon database.

####Alternative I/O Formats
Rather than passing RA (`-r`) and Dec (`-d`) separately, you can pass them to `--skypos` as follows.

    ./gFind.py -b 'NUV' --skypos '[176.919525856024,0.255696872807351]'

If you are interested only in the available raw exposure times and not the individual time ranges, pass the `--exponly` flag.

    ./gFind.py -b 'NUV' --skypos '[176.919525856024,0.255696872807351]' --exponly

###gAperture.py
_gAperture_ is the photometry tool.

####Lightcurve File Column Definitions
**[TBD]**

###gMap.py
_gMap_ is the image creation tool. It can generate integrated count, intensity, and response (equivalent to GALEX _cnt_, _int_ and _rrhr_) maps of arbitrary size<sup>+</sup>, shape and depth, including coadds across epochs and survey designation. It can also create "movie" (time-binned, multi-plane) versions of such maps.

_gMap_ writes all image files in the Flexible Image Transport System (FITS) standard, which is an uncompressed archival data format favored by many astronomy applications (and was used by the GALEX mission for archival products). FITS images generated by _gMap_ have headers which describe the World Coordinate System (WCS) information defining their orientation on the sky as well as effective exposure time and other metadata for the observation as a whole. For multi-frame images (movies), the per-plane information is described in a table in the FITS secondary HDU; this table describes start time, stop time, and effective exposure for each frame (within appropriately labelled columns).ß

<sup>+</sup>Caveat emptor. Arbitrariness is limited by your available RAM.

####Count Maps
Count (_cnt_) maps are integrated, aspect corrected but uncalibrated (that is, not adjusted for resonse or exposure time) images of the sky. Count images are good for "quick looks" at the data, to ensure that you are pointing in the location that you expected and that you are seeing the sources or features desired. But because they are not calibrated by either the relative response or the effective exposure time, you should not use them for photometric analyses of any kind.

You can create a count image from the command line by specifying the band (`-b`), sky position (`-r` and `-d` or just `--skypos`), the angular extent of the desired image in degrees (`--angle`), an output FITS filename (`--count`), and optionally a time range (`--t0` and `--t1` or just `--trange` or `--tranges`). Try the following simple command.

    ./gMap.py -b 'FUV' -r 176.919525856024 -d 0.255696872807351 --angle 0.5 --t0 766525332.995 --t1 766526576.995 --count 'count.fits'

Or, using the alternative formats for specifying sky position and time range, try the following.

    ./gMap.py -b 'FUV' --skypos '[176.919525856024,0.255696872807351]' --angle 0.5 --tranges '[[766525332.995, 766526576.995], [919755500.995, 919755916.995]]' --count 'count.fits'

Note that the most recent command created an image with two planes corresonding to the time ranges `[766525332.995, 766526576.995]` and `[919755500.995, 919755916.995]`, respectively. This behavior can be used to specify custom time ranges for a movie (see **Movies** below), or you can force _gMap_ to coadd these two time ranges with the `--coadd` flag.

    ./gMap.py -b 'FUV' --skypos '[176.919525856024,0.255696872807351]' --angle 0.5 --tranges '[[766525332.995, 766526576.995], [919755500.995, 919755916.995]]' --count 'count.fits' --coadd

However, if you do not specify any time range, _gMap_ will automatically use all available exposure. (And uses the same underlying function as _gFind_ to locate said exposure time, so therefore will also accept the `--gap` and `--minexp` and `--detsize` keywords if desired.)

    ./gMap.py -b 'FUV' --skypos '[176.919525856024,0.255696872807351]' --angle 0.5 --count 'count.fits' --coadd

####Response Maps
Relative response maps (_rrhr_) can be thought of as the detector flat field as projected onto the sky as a function of the detector boresight over time. Because the creation of relative response maps requires multiple computationally intensive interpolations, they take a long time to run. Intensity maps require response maps and therefore also take a long time to run. So, therefore, **WARNING: Making a response map will take longer than you expect. Do it sparingly.** If you want to create a response map, pass a FITS filename to the `--reponse` flag.

    ./gMap.py -b 'FUV' --skypos '[176.919525856024,0.255696872807351]' --angle 0.5 --response 'response.fits' --coadd

####Intensity Maps
Intensity maps (_int_) are integrated, aspect corrected images which have been corrected for both relative response as well as effective exposure time. If you need to perform photometric analysis on an image, you need an intensity map. Note the warning above about long run times for generating the _response_ for this. If you want to create an intensity map, pass a FITS filename to the `--intensity` flag.

    ./gMap.py -b 'FUV' --skypos '[176.919525856024,0.255696872807351]' --angle 0.5 --intensity 'intensity.fits' --coadd

####Movies
You can turn any of the above map types into a _movie_ (that is, a time-binned, multi-plane image), by passing the desired (raw) bin depth in seconds to the `--frame` parameter. For example, to create a count image with 100 second depth image planes, try the following.

    ./gMap.py -b 'FUV' -r 176.919525856024 -d 0.255696872807351 --angle 0.5 --t0 766525332.995 --t1 766526576.995 --count 'count.fits' --frame 100

##Common Questions, Issues, and Gotchas
1. **"My data is not available!"** You can verify that data for your desired target does or does not exist in the database and present by using the `gFind` commands described above. If data for your target is not available, there are two possible explanations: (1) we have not yet loaded those observations into the photon database, or (2) that target was never observed by the GALEX mission. As of this writing, we have only loaded about 5% of the total GALEX corpus into the photon database. You can confirm that your target was, indeed, observed by GALEX by searching for it in the [GALEX Catalog](http://galex.stsci.edu/GR6/?page=mastform). If your target was, indeed, observed by GALEX but has not yet been loaded into the gPhoton database, please contact Chase Million (chase.million@gmail.com) and Bernie Shiao (shiao@stsci.edu) with your target coordinates, and we will try to prioritize the associated data.
2. **Exposure Time.** When calculating the exposure time of a GALEX observation, one cannot simply subtract the end time from the start time. This _raw exposure time_ (t<sub>raw</sub>) must be adjusted into an _effective exposure time_ (t<sub>eff</sub>) which takes into account certain microchannel plate detector properties. In particular, the raw exposure time must be adjusted by both the _shutter_ (t<sub>shut</sub>) and the _deadtime ratio_ such that: t<sub>eff</sub> = (t<sub>raw</sub> - t<sub>shut</sub>) * deadtime.
    * **Shutter Correction.** GALEX did not observe all parts of the sky at all times. Even when GALEX was observing a particular part of the sky, there were times when little or no usable data was actually recorded. Data for specific targets or time ranges might exist but still be (practically) unusable for a variety for a variety of reasons which may include (1) that the spacecraft was not observing at nominal high voltage ("HVNOM"), (2) the spacecraft was observing a different part of the sky or in spectroscopic ("grism") mode, (3) a valid aspect solution could not be reconstructed from the available data (perhaps due to too few detectable reference stars), or (4) there was a temporary gap in the data due to a spacecraft or data transmissio anomaly. The _shutter correction_ accounts for these gaps by conceptualizing a "virtual shutter" that is considered closed whenever a gap of 0.05 seconds occurs in the data covering a particular source (as determined by boresight center and detector width) and within a particular time range. The integrated time of all such 0.05 second gaps in any given time range (t<sub>shut</sub>) is subtracted from t<sub>raw</sub>.
    * **Deadtime Correction** The GALEX microchannel plates could only detect a single event at a time; that is, when one photon event was being read into the electronics, any othe rincident photon events were missed entirely. This effectively reduced the exposure time of any observation by a ratio known as the "deadtime correction" which scaled as a function of field brightness (i.e. global count rate). Formally, the deadtime correction is the estimated fraction of the raw exposure time that the detector is not properly observing because it is engaged in readout. The deadtime changes from observation to observation and even moment to moment and, so, is ideally calculated continuously. There are two ways to estimate the deadtime.
        1. **Direct Measurement.**
        2. **Empirical Estimate.**
3. **Relative Response.** Microchannel plates do not have true "pixels" in the same sense as CCD or CMOS detectors. Microchannel plates are also prone to "burn in" or _local gain sag_ where regions of the detector exposured to high countrate sources can be local depleted of charge and suffer either temporary or permanent loss of sensitivity. To avoid this local gain sag, the telescope boresight did not remain stationary on the sky over the course of any observation, and this relatively motion was then accounted for in the _aspect correction_ stage of the calibration pipeline. **[...]**
4. **Gain Sag.** [COMING SOON]
