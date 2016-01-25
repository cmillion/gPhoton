"""
.. module:: CalibrationTools

   :synopsis: Something.
   @CHASE - elaborate on this please. Also, explain how this differs from
   CalUtils, which may be a source of confusion.@

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

# @CHASE - We need to avoid "import *", can we do specific imports here?@
import csv
import time
from astropy.io import fits as pyfits
import numpy as np
from FileUtils import *
import scipy.ndimage
from PhotonTools import *
from MCUtils import *
from CalUtils import *
from gnomonic import *
import cal

GPSSECS = 315532800+432000

# ------------------------------------------------------------------------------
def load_txy(csvfile):
    """
    Loads just the t,x,y columns from a photon CSV file.

    :param csvfile: Name of the photon event CSV file to read.

    :type csvfile: str

    :returns: tuple -- A four-element tuple containing arrays with the
    times, x position, y position, and flags from the photon event CSV file.
    """

    reader = csv.reader(open(csvfile, 'rb'), delimiter=',', quotechar='|')

    t, x, y, flags = [], [], [], []
    # Note, the columns are:
    # ['t', 'x', 'y', 'xa', 'ya', 'q', 'xi', 'eta', 'ra', 'dec', 'flags']

    for row in reader:
        t.append(float(row[0])/1000.)
        x.append(float(row[1]))
        y.append(float(row[2]))
        flags.append(float(row[10]))

    return np.array(t), np.array(x), np.array(y), np.array(flags, dtype='int16')
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def compute_deadtime(t, x, y, band, eclipse, trange=[[], []]):
    """
    Uses multiple methods to estimate the detector deadtime correction.
    The deadtime is an estimate of the fraction of time that the detector
    was unable to register new events because it was in the middle of
    readout. Deadtime should, therefore, not count as true exposure time.

    :param t: Set of photon event times.

    :type t: @CHASE - Is this a list, numpy.ndarray?@

    :param x: Set of photon event x positions.

    :type x: @CHASE - Is this a list, numpy.ndarray?@

    :param y: Set of photon event y positions.

    :type y: @CHASE - Is this a list, numpy.ndarray?@

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param eclipse: The eclipse number the data are taken from.

    :type eclipse: int @CHASE - Confirm this is int/long/float.@

    :param trange: The minimum and maximum times that define the range
    to calculate the dead time. If not supplied, the range is defined as the
    minimum and maximum times of the photon events.
    @CHASE - This is a list of lists, does it mean there can be multiple time
    ranges, if so, need to explain that in the description here, if not, then
    there should be some other default value used.@

    :type trange: list

    :returns: float or numpy.ndarray -- The dead time, as a percentage
    (0 < dead time < 1), during the entire observation or per time bin.
    @CHASE - please confirm return data type(s) and description.@
    """

    print "Computing deadtime correction..."
    refrate = 79.0 # counts per second
    minrate = refrate*.4
    maxrate = refrate+2.
    maxdiff = 0.1
    feeclkratio = 0.966
    tstep = 1. # seconds
    tec2fdead = 5.52e-6 # conversion from TEC to deadtime correction

    if not trange[0]:
        trange[0] = min(t)
    if not trange[1]:
        trange[1] = max(t)

    exptime = trange[1]-trange[0]

    stimt, stimx_as, stimy_as, stimix = find_stims(t, x, y, band, eclipse)
    print "Located "+str(len(stimt))+" stim events."

    dead = 0.
    # Compute using an emperical formula -- Method 1
    dead0 = tec2fdead*(len(t)/exptime)/feeclkratio
    print "	Simple correction w/ Method 1: "+str(dead0)

    # Toss more error checking into here esp. wrt minrate and maxrate
    if exptime <= tstep:
        tstep = 0.01
        bins = np.linspace(0., exptime-exptime%tstep, exptime//tstep+1)
        dead2 = (1.-((len(stimt)/exptime)/feeclkratio)/refrate)
    else:
        bins = np.linspace(0., exptime-exptime%tstep, exptime//tstep+1)
        h, xh = np.histogram(stimt-trange[0], bins=bins)
        ix = ((h <= maxrate) & (h >= minrate)).nonzero()[0]
        dead2 = (1.-((h[ix]/tstep)/feeclkratio)/refrate).mean()

    h, xh = np.histogram(t-trange[0], bins=bins)
    dead1 = (tec2fdead*(h/tstep)/feeclkratio).mean()
    print "	Correction w/ Method 1:        "+str(dead1)
    print "	Correction w/ Method 2:        "+str(dead2)

    # For short time slices, the "best" deadtime estimation method
    #  doesn't work very well.
    if exptime <= 5.:
        print "Short exposure. Using Method 1."
        dead = dead1
        if abs(dead1-dead0) > maxdiff:
            print "Warning: Deadtime corrections are suspect."
    else:
        print "Using Method 2."
        dead = dead2
        if abs(dead2-dead1) > maxdiff:
            print "Warning: Deadtime corrections are suspect."

    return dead
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def compute_shutter(t, trange=[[], []]):
    """
    Computes the detector shutter correction.
    The shutter correction accounts for short periods of time when no events
    were registered by the detector for any number of reasons. Any gap of
    longer than 0.05 seconds does not count as true exposure time.

    :param t: Set of photon event times.

    :type t: @CHASE - Is this a list, numpy.ndarray?@

    :param trange: The minimum and maximum times that define the range
    to calculate the shutter time. If not supplied, the range is defined as the
    minimum and maximum times of the photon events.
    @CHASE - This is a list of lists, does it mean there can be multiple time
    ranges, if so, need to explain that in the description here, if not, then
    there should be some other default value used.@

    :type trange: list

    :returns: float -- The total time lost due to shutter, i.e., the time spent
    while no events are registered by the detector.
    """

    if not trange[0]:
        trange[0] = min(t)
    if not trange[1]:
        trange[1] = max(t)

    exptime = trange[1]-trange[0]
    tstep = 0.05 # seconds

    bins = np.linspace(0., exptime-exptime%tstep, exptime//tstep+1)
    h, xh = np.histogram(t-trange[0], bins=bins)

    # If no counts are recorded for a tstep interval, consider the
    # virtual shutter to have been effectively closed during that time
    gaps = len((h == 0).nonzero()[0])

    return gaps*tstep
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def compute_exposure(t, x, y, flags, band, eclipse, trange=[[], []]):
    """
    Computes the effective, or 'true', exposure for the given data.

    :param t: Set of photon event times.

    :type t: @CHASE - Is this a list, numpy.ndarray?@

    :param x: Set of photon event x positions.

    :type x: @CHASE - Is this a list, numpy.ndarray?@

    :param y: Set of photon event y positions.

    :type y: @CHASE - Is this a list, numpy.ndarray?@

    :param flags: Set of flags for each photon event.

    :type flags: @CHASE - Is this a list, numpy.ndarray?@

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param eclipse: The eclipse number the data are taken from.

    :type eclipse: int @CHASE - Confirm this is int/long/float.@

    :param trange: The minimum and maximum times that define the range
    to calculate the effective exposure time. If not supplied, the range is
    defined as the minimum and maximum times of the photon events.
    @CHASE - This is a list of lists, does it mean there can be multiple time
    ranges, if so, need to explain that in the description here, if not, then
    there should be some other default value used.@

    :type trange: list

    :returns: float -- The effective exposure time, in seconds, during the
    specified time range. @CHASE - please confirm return data
    type(s) and description.@
    """

    # Use only unflagged data.
    # This should be done at the database level.
    ix = ((flags != 7) & (flags != 12)).nonzero()[0]

    if not len(ix):
        print "No unflagged data."
        return 0.
    if not trange[0]:
        trange[0] = min(t[ix])
    if not trange[1]:
        trange[1] = max(t[ix])

    exptime = trange[1]-trange[0]
    print "Gross exposure time is "+str(exptime)+" seconds."

    deadtime = exptime*compute_deadtime(t, x, y, band, eclipse,
                                        trange=trange)
    print "Removing "+str(deadtime)+" seconds of exposure from deadtime."

    ix = (flags == 0).nonzero()[0]
    shutter = compute_shutter(t[ix], trange=trange)

    if shutter:
        print ("Removing "+str(shutter)+" seconds of exposure from"
               " shutter.")
    print "Corrected exposure is "+str(exptime-deadtime-shutter)+" seconds."

    return exptime-deadtime-shutter
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def create_rr(csvfile, band, eclipse, aspfile=0., expstart=0., expend=0.,
              retries=20, detsize=1.25, pltscl=68.754932):
    """
    DEPRECATED: Creates a relative response map for an eclipse, given a
    photon list.

    :param csvfile: Name of CSV file containing photon events.

    :type csvfile: str

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param eclipse: The eclipse number the data are taken from.

    :type eclipse: int @CHASE - Confirm this is int/long/float.@

    :param aspfile: The aspect file to use with the CSV file.

    :type aspfile: str @CHASE - Please confirm this is a scalar string, and, if
    so, the default should be either '' or None instead of 0., with the
    checking of whether aspfile is specified updated accordingly.@

    :param expstart: Start of the exposure, in seconds.

    :type expstart: float

    :param expend: End of the exposure, in seconds.

    :type expend: float

    :param retries: Number of query retries before giving up.

    :type retries: int

    :param detsize: Effective size of the detector, in degrees.

    :type detsize: float

    :param pltscl: The plate scale in arcseconds.
    @CHASE - Confirm description and units are correct.@

    :type pltscl: float

    :returns: tuple -- A two-element tuple containing the relative response
    and the effective exposure time. The relative response is a 2D array of
    values.
    """

    aspum = pltscl/1000.0

    print "Loading flat file..."
    flat, flatinfo = cal.flat(band)
    npixx = flat.shape[0]
    npixy = flat.shape[1]
    pixsz = flatinfo['CDELT2']
    flatfill = detsize/(npixx*pixsz)

    print "Retrieving aspect data..."
    if aspfile:
        (aspra, aspdec, asptwist, asptime, aspheader,
         aspflags) = load_aspect([aspfile])
    else:
        (aspra, aspdec, asptwist, asptime, aspheader,
         aspflags) = web_query_aspect(eclipse, retries=retries)
    minasp = min(asptime)
    maxasp = max(asptime)
    print "			trange= ( "+str(minasp)+" , "+str(maxasp)+" )"
    ra0, dec0, roll0 = aspheader['RA'], aspheader['DEC'], aspheader['ROLL']
    print ("			[RA, DEC, ROLL] = ["+str(ra0)+", "+str(dec0)+
           ", "+str(roll0)+"]")

    print "Computing aspect vectors..."
    print "Calculating aspect solution vectors..."
    xi_vec, eta_vec = np.array([]), np.array([])
    xi_vec, eta_vec = gnomfwd_simple(ra0, dec0, aspra, aspdec, -asptwist,
                                     1.0/36000.0, 0.)

    # @CHASE - 'compute_flat_scale' seems to be undefined, need to import
    # 'galextools' in order to access this method, right?@
    flat_scale = compute_flat_scale(asptime.mean(), band)

    if not expstart:
        expstart = asptime.min()+GPSSECS
    if not expend:
        expend = asptime.max()+GPSSECS
    flatbuff = np.zeros([960, 960])
    # Rotate the flat into the correct orientation to start
    flatbuff[80:960-80, 80:960-80] = np.flipud(np.rot90(flat))
    expt = 0
    rr = np.zeros([960, 960])
    col = (((xi_vec/36000.)/(detsize/2.)*flatfill + 1.)/2. * npixx)-400.
    row = (((eta_vec/36000.)/(detsize/2.)*flatfill + 1.)/2. * npixy)-400.
    for i in xrange(len(asptime)-1):
        if (asptime[i]+GPSSECS) < expstart or (asptime[i]+GPSSECS) > expend:
            print "		", asptime[i]+GPSSECS, " out of range."
            continue
        elif (aspflags[i]%2 != 0) or (aspflags[i+1]%2 != 0):
            print "		", asptime[i]+GPSSECS, " flagged."
            continue
        else:
            rr += scipy.ndimage.interpolation.shift(
                scipy.ndimage.interpolation.rotate(flatbuff, -asptwist[i],
                                                   reshape=False, order=0,
                                                   prefilter=False), [col[i],
                                                                      row[i]],
                order=0, prefilter=False)
            expt += 1

    # Need to modify this to handle NUV files better.
    t, x, y, flags = load_txy(csvfile)
    exp = compute_exposure(t, x, y, flags, band, eclipse)
    deadt = compute_deadtime(t, x, y, band, eclipse)

    return rr*flat_scale*(1-deadt), exp
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def write_rr(csvfile, band, eclipse, rrfile, outfile, aspfile=0, expstart=0.,
             expend=0., exptime=0., imsz=960., retries=20):
    """
    Creates a relative response map for an eclipse, given a photon list
    file, and writes it to a FITS file.

    :param csvfile: Name of CSV file containing photon events.

    :type csvfile: str

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param eclipse: The eclipse number the data are taken from.

    :type eclipse: int @CHASE - Confirm this is int/long/float.@

    :param rrfile: The name of the existing relative response FITS file??
    @CHASE - I don't know what this is, since we are making a RR FITS file in
    this method, please update this description.@

    :type rrfile: str

    :param outfile: The name of the output relative response FITS file to make.

    :type outfile: str

    :param aspfile: The aspect file to use with the CSV file.

    :type aspfile: str @CHASE - Please confirm this is a scalar string, and, if
    so, the default should be either '' or None instead of 0., with the
    checking of whether aspfile is specified updated accordingly.@

    :param expstart: Start of the exposure, in seconds.

    :type expstart: float

    :param expend: End of the exposure, in seconds.

    :type expend: float

    :param exptime: @CHASE - This seems to be defined by the method 'create_rr',
    so can this actually be removed from the argument for this method?@

    :type exptime: float

    :param imsz: @CHASE - This doesn't seem to be used anywhere in the method,
    should remove (but need to make sure any calls elsewhere in the script don't
    supply this argument.@

    :type imsz: float

    :param retries: Number of query retries before giving up.

    :type retries: int
    """

    rr, exptime = create_rr(csvfile, band, eclipse, aspfile, expstart, expend,
                            retries=retries)

    if aspfile:
        (aspra, aspdec, asptwist, asptime, aspheader,
         aspflags) = load_aspect([aspfile])
    else:
        (aspra, aspdec, asptwist, asptime, aspheader,
         aspflags) = web_query_aspect(eclipse, retries=retries)

    ra0, dec0, roll0 = aspheader['RA'], aspheader['DEC'], aspheader['ROLL']

    if not expstart:
        mint = asptime.min()
    else:
        mint = expstart

    if not expend:
        maxt = asptime.max()
    else:
        maxt = expend

    hdulist = pyfits.open(rrfile)
    hdr = hdulist[0].header
    hdulist.close()
    hdr.update(key='expstart', value=mint)
    hdr.update(key='expend', value=maxt)
    if exptime:
        hdr.update(key='exptime', value=exptime)
    pyfits.writeto(outfile, rr, hdr, clobber=True)

    return
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def write_rrhr(rrfile, rrhrfile, outfile):
    """
    Turns a relative response (rr) into a high resolution relative response
    (rrhr) file with interpolation.

    :param rrfile: Relative Response FITS file to get header information from.
    @CHASE - Please correct/update this.@

    :type rrfile: str

    :param rrhrfile: High Resolution Relative Response FITS file to get header
    information from.
    @CHASE - Please correct/update this.@

    :type rrhrfile: str

    :param outfile: Name of the output FITS file to create.

    :type outfile: str
    """

    hdulist1 = pyfits.open(rrfile)
    hdr1 = hdulist1[0].header
    hdulist1.close

    hdulist0 = pyfits.open(rrhrfile)
    hdr0 = hdulist0[0].header
    hdulist0.close

    hdr0.update(key='expstart', value=hdr1['expstart'])
    hdr0.update(key='expend', value=hdr1['expend'])
    hdr0.update(key='exptime', value=hdr1['exptime'])

    rr = get_fits_data(rrfile)

    rrhr = scipy.ndimage.interpolation.zoom(rr, 4., order=0, prefilter=False)

    pyfits.writeto(outfile, rrhr, hdr0, clobber=True)

    return
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def write_int(cntfile, rrhrfile, oldint, outfile):
    """
    Writes out an intensity (int) map given a count (cnt) and a high
    resolution relative response (rrhr).

    :param cntfile: Name of count FITS file.
    @CHASE - Any more info that can be given here?@

    :type cntfile: str

    :param rrhrfile: Name of the High Resolution Relative Response FITS file.

    :type rrhrfile: str

    :param oldint: Name of Intensity file to get header from.
    @CHASE - Please correct/update this. Is this an existing file? Where did
    it come from if so?@

    :type oldint: str

    :param outfile: Name of the output FITS file to create.

    :type outfile: str
    """

    cnt = get_fits_data(cntfile)
    rrhr = get_fits_data(rrhrfile)

    hdulist1 = pyfits.open(cntfile)
    hdr1 = hdulist1[0].header
    hdulist1.close

    hdulist0 = pyfits.open(oldint)
    hdr0 = hdulist0[0].header
    hdulist0.close

    hdr0.update(key='expstart', value=hdr1['expstart'])
    hdr0.update(key='expend', value=hdr1['expend'])
    hdr0.update(key='exptime', value=hdr1['exptime'])

    int = cnt/rrhr
    int[np.where(np.isnan(int) == True)] = 0.

    pyfits.writeto(outfile, int, hdr0, clobber=True)

    return
# ------------------------------------------------------------------------------
