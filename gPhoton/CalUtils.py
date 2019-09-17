"""
.. module:: CalUtils
   :synopsis: Numerous methods for calibrating the raw photon event data.
       Many of these instantiate or make use of specific detector hardware
       parameters / constants related to the "static" or detector-space event
       calibration, including walk, wiggle, linearity, post-CSP, and stim
       scaling corrections.
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
from astropy.io import fits as pyfits
from builtins import str
from builtins import range
import csv
import numpy as np
import pandas as pd
# gPhoton imports.
import gPhoton.cal as cal
from gPhoton.galextools import isPostCSP
from gPhoton.MCUtils import rms, print_inline, get_fits_header, get_tbl_data

# ------------------------------------------------------------------------------
def clk_cen_scl_slp(band, eclipse):
    """
    Return the detector clock, center, scale, and slope constants. These are
        empirically determined detector-space calibration parameters that help
        define the relationship between raw event data and its physical position
        on the detector.

    :param band: The band to return the constants for, either 'FUV' or 'NUV'.

    :type band: str

    :param eclipse: The eclipse number to return the constants for.

    :type eclipse: int

    :returns: tuple -- a tuple containing the x-clock, y-clock, x-center,
        y-center, x-scale, y-scale, x-slope, and y-slope constants.
    """

    band = band.upper()

    if band == 'FUV':
        xclk, yclk = 1997., 1993.
        xcen, ycen = 7200., 6670.
        xscl, yscl = 7.78, 10.73
        xslp, yslp = 0., 0.
    elif band == 'NUV':
        xclk, yclk = 2007., 1992.
        # Special value for post-CSP event.
        if eclipse >= 38150:
            yclk = 2016.
        xcen, ycen = 7400., 6070.
        xscl, yscl = 8.79, 14.01
        xslp, yslp = 0.53, 0.
    else:
        # Might need to raise a ValueError or some other proper error?
        print('Band must be either fuv or nuv ... Exiting.')
        return

    return xclk, yclk, xcen, ycen, xscl, yscl, xslp, yslp
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def xieta2colrow(xi, eta, detsize, fill, npixx, npixy):
    """
    Convert detector xi and eta values to image column and row.

    :param xi: One coordinate of the unprojected detector position.

    :type xi: numpy.ndarray

    :param eta: Second coordinate of the unprojected detector position.

    :type eta: numpy.ndarray

    :param fill: Ration of the detector extent to the image extent.

    :type fill: float

    :param npixx: Number of pixels in the x dimension of an image.

    :type npixx: int

    :param nipxy: Number of pixels in the y dimension of an image.

    :type npixy: int

    :returns: tuple -- A tuple containing the image column and row values.
    """

    col = (((xi/36000.)/(detsize/2.)*fill+1.)/2.*npixx)

    row = (((eta/36000.)/(detsize/2.)*fill+1.)/2.*npixy)

    return col, row
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def avg_stimpos(band, eclipse):
    """
    Define the mean detector stim positions.

    :param band: The band to return the average stim positions for,
        either 'FUV' or 'NUV'.

    :type band: str

    :param eclipse: The eclipse number to return the average stim positions for.

    :type eclipse: int

    :returns: dict -- A dict containing the average x and y positions of the
        four stims (for the requested band).
    """

    if band == 'FUV':
        avgstim = ({'x1':-2541.88, 'x2':2632.06, 'x3':-2541.53,
                    'x4':2631.68, 'y1':2455.28, 'y2':2455.02,
                    'y3':-2550.89, 'y4':-2550.92})

    elif band == 'NUV':
        if eclipse >= 38268:
            # The average stim positions after the clock change (post-CSP).
            avgstim = ({'x1':-2722.53, 'x2':2470.29, 'x3':-2721.98,
                        'x4':2471.09, 'y1':2549.96, 'y2':2550.10,
                        'y3':-2538.57, 'y4':-2538.62})
        else:
            # The average stim positions for pre-CSP data (eclipse 37423).
            avgstim = ({'x1':-2722.27, 'x2':2468.84, 'x3':-2721.87,
                        'x4':2469.85, 'y1':2453.89, 'y2':2453.78,
                        'y3':-2565.81, 'y4':-2567.83})

    else:
        # Might need to raise a ValueError or other proper error here.
        print("Error: No valid band specified.")

    return avgstim
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def find_stims_index(x, y, band, eclipse, margin=90.001):
    """
    Given a list of detector x,y positions of events, returns four
        arrays that contain the indices of likely stim events for that stim,
        i.e., the first array contains positions for stim1, the second array has
        positions of stim2, etc.

        Example of how the return indexes are used: x[index1], y[index1] would
        give all of the event positions for stim1.

    :param x: Detector 'x' positions to identify likely stim events from.

    :type x: list

    :param y: Detector 'y' positions to identify likely stim events from.

    :type y: list

    :param band: The band to return the constants for, either 'FUV' or 'NUV'.

    :type band: str

    :param eclipse: The eclipse number to return the average stim positions for.

    :type eclipse: int

    :param margin: +/- extent in arcseconds defining search box

    :type margin: float

    :returns: tuple -- A four-element tuple containing arrays of indexes that
        correpond to the event positions for stim1, stim2, stim3, and stim4.
    """

    # [Future] This method could be programmed better. Consider using numpy
    # "where" and the logical '&' operator, instead of .nonzero()?

    # Plate scale (in arcsec/mm).
    pltscl = 68.754932
    # Plate scale in arcsec per micron.
    aspum = pltscl/1000.0
    # [Future] Could define the plate scale (in both units) as a constant for
    # the module, since it is used in many places?

    x_as = np.array(x)*aspum
    y_as = np.array(y)*aspum

    avg = avg_stimpos(band, eclipse)

    index1 = (
        ((x_as > (avg['x1']-margin)) & (x_as < (avg['x1']+margin)) &
         (y_as > (avg['y1']-margin)) & (y_as < (avg['y1']+margin))).nonzero()[0]
        )
    index2 = (
        ((x_as > (avg['x2']-margin)) & (x_as < (avg['x2']+margin)) &
         (y_as > (avg['y2']-margin)) & (y_as < (avg['y2']+margin))).nonzero()[0]
        )
    index3 = (
        ((x_as > (avg['x3']-margin)) & (x_as < (avg['x3']+margin)) &
         (y_as > (avg['y3']-margin)) & (y_as < (avg['y3']+margin))).nonzero()[0]
        )
    index4 = (
        ((x_as > (avg['x4']-margin)) & (x_as < (avg['x4']+margin)) &
         (y_as > (avg['y4']-margin)) & (y_as < (avg['y4']+margin))).nonzero()[0]
        )

    return index1, index2, index3, index4
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def stimcount(data, band, t0, t1, margin=90.001):
    """
    Given a dict() that contains a list of 't' event times and 'x' and 'y'
        detector positions, returns the total number of stim events within a
        time range.

    :param data: The event times and detector positions from which to count the
        total number of stim events.

    :type data: dict

    :param band: The band to return the constants for, either 'FUV' or 'NUV'.

    :type band: str

    :param t0: The mininum of the time range to report stim counts within.

    :type t0: float

    :param t1: The maximum of the time range to report stim counts within.

    :type t1: float

    :param margin: +/- extent in arcseconds defining search box

    :type margin: float

    :returns: int -- The total number of stim events.
    """

    # [Future] This method duplicates major functionality of the "" method.
    # Consider sharing duplicate calculations (particularly finding stim event
    # locations) as a shared method?

    # Plate scale (in arcsec/mm).
    pltscl = 68.754932
    # Plate scale in arcsec per micron.
    aspum = pltscl/1000.0
    # [Future] If these are made module constants, can remove from this method.

    eclipse = 100. if isPostCSP(t0) else 40000. # HACK: for backwards comp.

    avgstim = avg_stimpos(band, eclipse)

    time = ((np.array(data['t']) >= t0) & (np.array(data['t']) <= t1))

    stim1 = ((np.array(data['x']) > (avgstim['x1']-margin)) &
             (np.array(data['x']) < (avgstim['x1']+margin)) &
             (np.array(data['y']) > (avgstim['y1']-margin)) &
             (np.array(data['y']) < (avgstim['y1']+margin)))
    stim2 = ((np.array(data['x']) > (avgstim['x2']-margin)) &
             (np.array(data['x']) < (avgstim['x2']+margin)) &
             (np.array(data['y']) > (avgstim['y2']-margin)) &
             (np.array(data['y']) < (avgstim['y2']+margin)))
    stim3 = ((np.array(data['x']) > (avgstim['x3']-margin)) &
             (np.array(data['x']) < (avgstim['x3']+margin)) &
             (np.array(data['y']) > (avgstim['y3']-margin)) &
             (np.array(data['y']) < (avgstim['y3']+margin)))
    stim4 = ((np.array(data['x']) > (avgstim['x4']-margin)) &
             (np.array(data['x']) < (avgstim['x4']+margin)) &
             (np.array(data['y']) > (avgstim['y4']-margin)) &
             (np.array(data['y']) < (avgstim['y4']+margin)))

    ix = np.where(time & (stim1 | stim2 | stim3 | stim4))

    return len(ix[0])
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def totalcount(data, t0, t1):
    """
    Given a dict() containg 't', a list of global even times, return
        the total number of events within the time range.

    :param data: The event times and detector positions from which to count the
        total number of photon events.

    :type data: dict

    :param t0: The mininum of the time range to report total counts within.

    :type t0: float

    :param t1: The maximum of the time range to report total counts within.

    :type t1: float

    :returns: int -- The total number of events within the given time range.
    """
    ix = np.where((np.array(data['t']) >= t0) & (np.array(data['t']) <= t1))

    return len(ix[0])
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def deadtime_method0(data, t0, t1, band, feeclkratio=0.966, tec2fdead=5.52e-6):
    """
    Given a dict() containing 't', a list of global event times, computes
        the deadtime using an empirical formula based on global count rate over
        the whole time range.

    :param data: The event times from which to calculate the deadtime.

    :type data: dict

    :param t0: The mininum of the time range to calculate the deadtime within.

    :type t0: float

    :param t1: The maximum of the time range to calculate the deadtime within.

    :type t1: float

    :param band: The band to return the constants for, either 'FUV' or 'NUV'.
        Not actually used in this function, but retained to make this function
        definition consistent with deadtime_method[12].

    :type band: str

    :param feeclkratio: Ratio of Front End Electronics clock rates.

    :type feeclkratio: float

    :param tec2fdead: The nominal amount of time following an event that the
        detector is unable to detect another event.

    :type tec2fdead: float

    :returns: float -- The deadtime for the specified time range.
    """
    exptime = t1-t0

    totcount = totalcount(data, t0, t1)

    return tec2fdead*(totcount/exptime)/feeclkratio
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def deadtime_method1(data, t0, t1, band, feeclkratio=0.966, tec2fdead=5.52e-6,
                     tstep=1.):
    """
    Given a dict() containing 't', a list of global event times, computes
        the deadtime using an empirical formula based on global count rates, put
        into bins of depth equal to `tstep` seconds and averaged.

    :param data: The event times and detector positions from which to calculate
        the deadtime.

    :type data: dict

    :param t0: The mininum of the time range to calculate the deadtime within.

    :type t0: float

    :param t1: The maximum of the time range to calculate the deadtime within.

    :type t1: float

    :param band: The band to return the constants for, either 'FUV' or 'NUV'.
        Not actually used in this function, but retained to make this function
        definition consistent with deadtime_method[02].

    :type band: str

    :param feeclkratio: Ratio of Front End Electronics clock rates.

    :type feeclkratio: float

    :param tec2fdead: The nominal amount of time following an event that the
        detector is unable to detect another event.

    :type tec2fdead: float

    :type tec2fdead: float

    :param tstep: Bin size (in seconds) to use when reporting the deadtime.

    :type tstep: float

    :returns: numpy.ndarray -- The deadtimes within the specified time range,
        split into the requested bin size.
    """

    exptime = t1-t0

    bins = np.linspace(0., exptime-exptime%tstep, exptime//tstep+1)+t0

    h = np.zeros(len(bins))

    for i, t in enumerate(bins):
        h[i] = totalcount(data, t, t+tstep-0.0001)

    return (tec2fdead*(h/tstep)/feeclkratio).mean()
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def deadtime_method2(data, t0, t1, band, refrate=79., feeclkratio=0.966,
                     tstep=1., refrange=[.4, 2.]):
    """
    Given a list of global event times, computes the deadtime through
        direct comparison of the stim rate to the reference rate in bins of
        depth `tstep` seconds, and trimmed of outliers. This is close to the
        deadtime method used by the mission pipeline.

    :param data: The event times and 'x' and 'y' detector positions from
        which to calculate the deadtime.

    :type data: dict

    :param t0: The mininum of the time range to calculate the deadtime within.

    :type t0: float

    :param t1: The maximum of the time range to calculate the deadtime within.

    :type t1: float

    :param band: The band to return the constants for, either 'FUV' or 'NUV'.

    :type band: str

    :param refrate: The reference rate of stims in counts/sec.

    :type refrate: float

    :param feeclkratio: Ratio of Front End Electronics clock rates.

    :type feeclkratio: float

    :param tstep: Bin size (in seconds) to use when reporting the deadtime.

    :type tstep: float

    :param refrange: minimum and maximum multiplicative of the reference stim
        rate that will be considered a legitimate / valid measurement

    :type refrange: list

    :returns: numpy.ndarray -- The deadtimes within the specified time range,
        split into the requested bin size.
    """

    eclipse = 100. if isPostCSP(t0) else 40000. # HACK: for backwards comp.

    exptime = t1-t0

    bins = np.linspace(0., exptime-exptime%tstep, exptime//tstep+1)+t0

    h = np.zeros(len(bins))

    for i, t in enumerate(bins):
        h[i] = stimcount(data, band, t, t+tstep-0.0001, eclipse)/exptime

    (minrate, maxrate) = (refrate*refrange[0], refrate+refrange[1])

    ix = ((h <= maxrate) & (h >= minrate)).nonzero()[0]

    return (1.-((h[ix]/tstep)/feeclkratio)/refrate).mean()
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def deadtime_fromlist(photonfile, t0, t1, band, method=0,
                      colnames=['t', 'x', 'y', 'xa', 'ya', 'q', 'xi', 'eta',
                                'ra', 'dec', 'flags']):
    """
    Given a photon event list, calculate the deadtime using the specified
        method.

    :param photonfile: Name of the file containing photon events.

    :type photonfile: str

    :param t0: The mininum of the time range to calculate the deadtime within.

    :type t0: float

    :param t1: The maximum of the time range to calculate the deadtime within.

    :type t1: float

    :param band: The band to return the constants for, either 'FUV' or 'NUV'.

    :type band: str

    :param method: which deadtime_method[012] to use when calculating deadtime

    :type method: int

    :param colnames: The column names within the photon event file.

    :type colnames: list

    :returns: float or numpy.ndarray -- The deadtime within the requested time
        range. If using Method 0, return value is a float, otherwise, return
        value is an array containing the deadtimes per bin.
    """

    # Read photon event file.
    data = pd.io.parsers.read_csv(photonfile, names=colnames)

    # Calculate deadtime using the requested method.
    return {0:deadtime_method0,
            1:deadtime_method1,
            2:deadtime_method2}[method](data, t0, t1, band)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def find_stims(t, x, y, band, eclipse):
    """
    Returns t,x,y and identity (1-4) of likely stims in the input arrays.

    :param t: Photon event times within which to search for stim events.

    :type t: numpy.ndarray

    :param x: Photon event detector x positions within which to search for stim
        events.

    :type x: numpy.ndarray

    :param y: Photon event detector y positions within which to search for stim
        events.

    :type y: numpy.ndarray

    :param band: The band to return the constants for, either 'FUV' or 'NUV'.

    :type band: str

    :param eclipse: The eclipse number to return the constants for.

    :type eclipse: int

    :returns: tuple -- A four-element tuple containing arrays of the times,
        detector x, detector y, and stim ID (1, 2, 3 or 4) of the likely stim
        events.
    """

    # [Future] This method could be programmed better. Consider using numpy
    # "where" and the logical '&' operator, instead of .nonzero()?

    # Plate scale (in arcsec/mm).
    pltscl = 68.754932
    # Plate scale in arcsec per micron.
    aspum = pltscl/1000.0
    # [Future] If these are made module constants, can remove from this method.

    stimt, stimix = np.array([]), np.array([])
    stimx_as, stimy_as = np.array([]), np.array([])

    x_as = np.array(x)*aspum
    y_as = np.array(y)*aspum

    avg = avg_stimpos(band, eclipse)

    index = (
        ((x_as > (avg['x1']-90.001)) & (x_as < (avg['x1']+90.001)) &
         (y_as > (avg['y1']-90.001)) & (y_as < (avg['y1']+90.001))).nonzero()[0]
        )

    stimt = np.append(stimt, t[index])
    stimx_as = np.append(stimx_as, x_as[index])
    stimy_as = np.append(stimy_as, y_as[index])
    stimix = np.append(stimix, np.zeros(len(index))+1)

    index = (
        ((x_as > (avg['x2']-90.001)) & (x_as < (avg['x2']+90.001)) &
         (y_as > (avg['y2']-90.001)) & (y_as < (avg['y2']+90.001))).nonzero()[0]
        )

    stimt = np.append(stimt, t[index])
    stimx_as = np.append(stimx_as, x_as[index])
    stimy_as = np.append(stimy_as, y_as[index])
    stimix = np.append(stimix, np.zeros(len(index))+2)

    index = (
        ((x_as > (avg['x3']-90.001)) & (x_as < (avg['x3']+90.001)) &
         (y_as > (avg['y3']-90.001)) & (y_as < (avg['y3']+90.001))).nonzero()[0]
        )

    stimt = np.append(stimt, t[index])
    stimx_as = np.append(stimx_as, x_as[index])
    stimy_as = np.append(stimy_as, y_as[index])
    stimix = np.append(stimix, np.zeros(len(index))+3)

    index = (
        ((x_as > (avg['x4']-90.001)) & (x_as < (avg['x4']+90.001)) &
         (y_as > (avg['y4']-90.001)) & (y_as < (avg['y4']+90.001))).nonzero()[0]
        )

    stimt = np.append(stimt, t[index])
    stimx_as = np.append(stimx_as, x_as[index])
    stimy_as = np.append(stimy_as, y_as[index])
    stimix = np.append(stimix, np.zeros(len(index))+4)

    index = np.argsort(stimt)

    return stimt[index], stimx_as[index], stimy_as[index], stimix[index]
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def get_stim_coefs(ssdfile):
    """
    Computes the stim scaling coefficients based on a ssdfile (which contains
        information on the on-detector spatial separation of stims as a function
        of time).

    :param ssdfile: The name of the Stim Separatation Data (SSD) file.

    :type ssdfile: str

    :returns: tuple -- A two-element tuple containing the stim scaling
        coefficients.
    """
    tbl = get_tbl_data(ssdfile)

    c11 = sum(tbl[:, 2])
    c12 = sum(tbl[:, 0]*tbl[:, 2])
    c13 = sum(tbl[:, 1]*tbl[:, 2])
    c22 = sum(tbl[:, 0]*tbl[:, 0]*tbl[:, 2])
    c23 = sum(tbl[:, 1]*tbl[:, 2]*tbl[:, 0])

    stim_coef1 = ((c13*c12)-(c23*c11))/((c12*c12)-(c22*c11))
    stim_coef0 = (c13-(c12*stim_coef1))/c11

    return stim_coef0, stim_coef1
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def find_fuv_offset(scstfile):
    """
    Computes NUV->FUV center offset based on a lookup table.

    :param scstfile: Name of the spacecraft state (-scst) FITS file.

    :type scstfile: str

    :returns: tuple -- A two-element tuple containing the x and y offsets.
    """

    fodx_coef_0, fody_coef_0, fodx_coef_1, fody_coef_1 = (0., 0., 0., 0.)

    scsthead = get_fits_header(scstfile)

    print("Reading header values from scst file: ", scstfile)

    try:
        eclipse = int(scsthead['eclipse'])
    except:
        print("WARNING: ECLIPSE is not defined in SCST header.")
        try:
            eclipse = int(photon_file.split('/')[-1].split('-')[0][1:])
            print("         Using {e} from filename.".format(e=eclipse))
        except:
            print("         Unable to infer eclipse from filename.")
            return 0., 0.

    try:
        fdttdc = float(scsthead['FDTTDC'])
    except KeyError:
        print("WARNING: FUV temperature value missing from SCST.")
        print("         This is probably not a valid FUV observation.")
        #raise
        return 0., 0.

    print("Offsetting FUV image for eclipse {e} at {t} degrees.".format(
        e=eclipse, t=fdttdc))

    fodx_coef_0 = cal.offset('x')[eclipse-1, 1]
    fody_coef_0 = cal.offset('y')[eclipse-1, 1]

    fodx_coef_1 = 0.
    fody_coef_1 = 0.3597

    if (fdttdc <= 20.) or (fdttdc >= 40.):
        print("ERROR: FDTTDC is out of range at {t}".format(t=fdttdc))
        return 0., 0.
    else:
        xoffset = fodx_coef_0 - (fodx_coef_1 * (fdttdc - 29.))
        yoffset = fody_coef_0 - (fody_coef_1 * (fdttdc - 29.))
        print("Setting FUV offsets to x={x}, y={y}".format(x=xoffset,
                                                           y=yoffset))

    return xoffset, yoffset
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def post_csp_caldata():
    """
    Loads the calibration data for after the CSP event.

    :returns: tuple -- A six-element tuple containing the wiggle, walk, and
        clock corrections. See the calibration paper for details.
    """

    print("Loading post-CSP wiggle file...")

    wig2fits, wig2head = cal.wiggle2()
    wig2 = np.zeros([128, 8, 32, 8])

    for i in range(len(wig2fits)):
        ya = wig2fits[i][0]
        yb = wig2fits[i][1]
        xb = wig2fits[i][2]
        yy = wig2fits[i][3]
        ycor = wig2fits[i][4]
        wig2[yy][yb][ya][xb] = ycor

    wig2data = {'start':wig2head['Y_AS_0'], 'inc':wig2head['Y_AS_INC']}

    print("Loading post-CSP walk file...")

    wlk2fits, wlk2head = cal.walk2()
    wlk2 = np.zeros([100, 8, 32])

    for i in range(len(wlk2fits)):
        q = wlk2fits[i][0]
        yb = wlk2fits[i][1]
        yy = wlk2fits[i][2]
        ycor = wlk2fits[i][3]
        wlk2[yy][yb][q] = ycor

    wlk2data = {'start':wlk2head['Y_AS_0'], 'inc':wlk2head['Y_AS_INC']}

    print("Loading post-CSP clock file...")

    clk2fits, clk2head = cal.clock2()
    clk2 = np.zeros([100, 8])

    for i in range(len(clk2fits)):
        yb = clk2fits[i][0]
        yy = clk2fits[i][1]
        ycor = clk2fits[i][2]
        clk2[yy][yb] = ycor

    clk2data = {'start':clk2head['Y_AS_0'], 'inc':clk2head['Y_AS_INC']}

    return wig2, wig2data, wlk2, wlk2data, clk2, clk2data
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def rtaph_yap(ya, yb, yamc):
    """
    For post-CSP data, 'wrap' the YA value for YA in [0,1]. From rtaph.c.

    :param ya: Y axis wiggle.

    :type ya: numpy.ndarray

    :param yb: Y axis coarse clock.

    :type yb: numpy.ndarray

    :param yamc: Raw Y detector position in FEE pixels.

    :type yamc: numpy.ndarray

    :returns: numpy.ndarray
    """

    yap = np.append([], ya)
    ix = ((yb > 1) & (yb < 5)).nonzero()[0]

    ix1 = ((ya[ix] == 0) & (yamc[ix] > -50)).nonzero()[0]
    yap[ix[ix1]] += 32

    ix1 = ((ya[ix] == 1) & (yamc[ix] > -10)).nonzero()[0]
    yap[ix[ix1]] += 32

    yap = np.array(yap, dtype='int64') % 32

    return yap
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def rtaph_yac(yactbl, ya, yb, yamc, eclipse):
    """
    Compute a the YAmC (yac) correction to the Y detector FEE position.

    :param yactbl: yac correction lookup table

    :type yactbl: numpy.ndarry

    :param ya: Y axis wiggle.

    :type ya: numpy.ndarray

    :param yb: Y axis coarse clock.

    :type yb: numpy.ndarray

    :param yamc: Raw Y detector position in FEE pixels.

    :type yamc: numpy.ndarray

    :param eclipse: The eclipse number to return the constants for.

    :type eclipse: int

    :returns: numpy.ndarray
    """

    yac = np.zeros(len(ya))
    if eclipse <= 37460:
        return yac

    yap = rtaph_yap(ya, yb, yamc)

    yac = yactbl[yap, yb]

    return yac
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def rtaph_yac2(q, xb, yb, ya, y, aspum, wig2, wig2data, wlk2, wlk2data, clk2,
               clk2data):
    """
    Compute a secondary correction to the YAmC (yac) detector FEE Y position.

    :param q: Detected pulse height

    :type q: numpy.ndarray

    :param xb: X axis coarse clock.

    :type xb: numpy.ndarray

    :param yb: Y axis coarse clock.

    :type yb: int

    :param ya: Y axis wiggle.

    :type ya: int

    :param y: Detector y position.

    :type y: numpy.ndarray

    :param aspum: Detector arcseconds per micrometer.

    :type aspum: float

    :param wig2: Secondary wiggle correction lookup table.

    :type wig2: numpy.ndarray

    :param wig2data: Secondary wiggle reference values.

    :type wig2data: dict

    :param wlk2: Secondary walk correction lookup table.

    :type wlk2: numpy.ndarray

    :param wlk2data: Secondary walk correction reference values.

    :type wlk2data: dict

    :param clk2: Secondary clock correction lookup table.

    :type clk2: numpy.ndarray

    :param clk2data: Secondary clock correction reference values.

    :type clk2data: dict

    :returns: numpy.ndarray -- Secondary YAmC corrections.
    """

    yac = 0
    y_as = y*aspum
    yac_as = np.zeros(len(y_as))
    ix = ((y_as > -2000)&(y_as < 2000)).nonzero()[0]

    ii = (np.array(y_as, dtype='int64')-wig2data['start'])/wig2data['inc']
    yac_as[ix] = wig2[np.array(ii[ix],dtype='int64'), yb[ix], ya[ix], xb[ix]]

    ii = (np.array(y_as, dtype='int64')-wlk2data['start'])/wlk2data['inc']
    yac_as[ix] = wlk2[np.array(ii[ix],dtype='int64'), yb[ix], q[ix]]

    ii = (np.array(y_as, dtype='int64')-clk2data['start'])/clk2data['inc']
    yac_as[ix] = clk2[np.array(ii[ix],dtype='int64'), yb[ix]]

    return yac_as/aspum
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def raw6_to_stims(raw6file, band, eclipse, margin=90.001):
    """
    Extracts stim events from a raw6 file.

    :param raw6file: The name of the raw6 FITS file to read.

    :type raw6file: str

    :param band: The band to return the stim data for, either 'FUV' or 'NUV'.

    :type band: str

    :param eclipse: The eclipse number to return the stim data for.

    :type eclipse: int

    :param margin: +/- extent in arcseconds defining stim search box

    :type margin: float

    :returns: tuple -- A four-element tuple containing data from each stim. The
        data from each stim are stored in dicts.
    """

    print("Extracting stim data from ", raw6file, " ...")
    print("         Using a search box with sides of ", margin, " arcseconds.")

    # This is unscoped for some reason... so I'm just coding it.
    (xclk, yclk, xcen, ycen, xscl, yscl, xslp, yslp) = clk_cen_scl_slp(band,
                                                                       eclipse)

    chunksz = 1000000
    print("Loading raw6 file...")
    raw6hdulist = pyfits.open(raw6file, memmap=1)
    raw6htab = raw6hdulist[1].header
    nphots = raw6htab['NAXIS2']

    stim1 = (
        {'t':np.array([]), 'q':np.array([]), 'xb':np.array([]),
         'xamc':np.array([]), 'yamc':np.array([]), 'xa':np.array([]),
         'ya':np.array([]), 'x':np.array([]), 'y':np.array([]),
         'yb':np.array([]), 'yap':np.array([])}
        )
    stim2 = (
        {'t':np.array([]), 'q':np.array([]), 'xb':np.array([]),
         'xamc':np.array([]), 'yamc':np.array([]), 'xa':np.array([]),
         'ya':np.array([]), 'x':np.array([]), 'y':np.array([]),
         'yb':np.array([]), 'yap':np.array([])}
        )
    stim3 = (
        {'t':np.array([]), 'q':np.array([]), 'xb':np.array([]),
         'xamc':np.array([]), 'yamc':np.array([]), 'xa':np.array([]),
         'ya':np.array([]), 'x':np.array([]), 'y':np.array([]),
         'yb':np.array([]), 'yap':np.array([])}
        )
    stim4 = (
        {'t':np.array([]), 'q':np.array([]), 'xb':np.array([]),
         'xamc':np.array([]), 'yamc':np.array([]), 'xa':np.array([]),
         'ya':np.array([]), 'x':np.array([]), 'y':np.array([]),
         'yb':np.array([]), 'yap':np.array([])}
        )

    print("")

    for i in range(int(nphots/chunksz)+1):
        csvrows = []
        chunkbeg, chunkend = i*chunksz, (i+1)*chunksz-1
        if chunkend > nphots:
            chunkend = nphots-1
        chunkid = " "+str(i+1)+" of "+str(int(nphots/chunksz)+1)+": "
        print_inline(chunkid+"Unpacking raw6 data...")
        t = np.array(raw6hdulist[1].data.field('t')[chunkbeg:chunkend])
        phb1 = np.array(
            raw6hdulist[1].data.field('phb1')[chunkbeg:chunkend], dtype='int64')
        phb2 = np.array(
            raw6hdulist[1].data.field('phb2')[chunkbeg:chunkend], dtype='int64')
        phb3 = np.array(
            raw6hdulist[1].data.field('phb3')[chunkbeg:chunkend], dtype='int64')
        phb4 = np.array(
            raw6hdulist[1].data.field('phb4')[chunkbeg:chunkend], dtype='int64')
        phb5 = np.array(
            raw6hdulist[1].data.field('phb5')[chunkbeg:chunkend], dtype='int64')

        q = ((phb4 & 3) << 3) + ((phb5 & 224) >> 5)
        xb = phb1 >> 5
        xamc = (np.array(((phb1 & 31) << 7), dtype='int16') +
                np.array(((phb2 & 254) >> 1), dtype='int16') -
                np.array(((phb1 & 16) << 8), dtype='int16'))
        yb = ((phb2 & 1) << 2) + ((phb3 & 192) >> 6)
        yamc = (np.array(((phb3 & 63) << 6), dtype='int16') +
                np.array(((phb4 & 252) >> 2), dtype='int16') -
                np.array(((phb3 & 32) << 7), dtype='int16'))
        xa = ((phb5 & 16) >> 4) + ((phb5 & 3) << 3) + ((phb5 & 12) >> 1)
        xraw0 = xb*xclk + xamc
        yraw0 = yb*yclk + yamc
        ya = (np.array(((((yraw0/(2*yclk) - xraw0/(2*xclk)) + 10)*32) + xa),
                       dtype='int64') % 32)
        xraw = (xraw0 + np.array((((xa+7) % 32) - 16), dtype='int64') * xslp)
        yraw = (yraw0 + np.array((((ya+7) % 32) - 16), dtype='int64') * yslp)
        x = (xraw - xcen)*xscl
        y = (yraw - ycen)*yscl

        index1, index2, index3, index4 = find_stims_index(x, y, band, eclipse,
                                                          margin)

        # [Future] There may well be a better way to do these assignments.
        stim1['t'] = np.append(stim1['t'], t[index1])
        stim1['x'] = np.append(stim1['x'], x[index1])
        stim1['y'] = np.append(stim1['y'], y[index1])
        stim1['q'] = np.append(stim1['q'], q[index1])
        stim1['xa'] = np.append(stim1['xa'], xa[index1])
        stim1['xb'] = np.append(stim1['xb'], ya[index1])
        stim1['ya'] = np.append(stim1['ya'], ya[index1])
        stim1['yb'] = np.append(stim1['yb'], yb[index1])
        stim1['xamc'] = np.append(stim1['xamc'], xamc[index1])
        stim1['yamc'] = np.append(stim1['yamc'], yamc[index1])
        stim1['yap'] = np.append(stim1['yap'], rtaph_yap(ya[index1],
                                                         yb[index1],
                                                         yamc[index1]))
        stim2['t'] = np.append(stim2['t'], t[index2])
        stim2['x'] = np.append(stim2['x'], x[index2])
        stim2['y'] = np.append(stim2['y'], y[index2])
        stim2['q'] = np.append(stim2['q'], q[index2])
        stim2['xa'] = np.append(stim2['xa'], xa[index2])
        stim2['xb'] = np.append(stim2['xb'], ya[index2])
        stim2['ya'] = np.append(stim2['ya'], ya[index2])
        stim2['yb'] = np.append(stim2['yb'], yb[index2])
        stim2['xamc'] = np.append(stim2['xamc'], xamc[index2])
        stim2['yamc'] = np.append(stim2['yamc'], yamc[index2])
        stim2['yap'] = np.append(stim2['yap'], rtaph_yap(ya[index2],
                                                         yb[index2],
                                                         yamc[index2]))
        stim3['t'] = np.append(stim3['t'], t[index3])
        stim3['x'] = np.append(stim3['x'], x[index3])
        stim3['y'] = np.append(stim3['y'], y[index3])
        stim3['q'] = np.append(stim3['q'], q[index3])
        stim3['xa'] = np.append(stim3['xa'], xa[index3])
        stim3['xb'] = np.append(stim3['xb'], ya[index3])
        stim3['ya'] = np.append(stim3['ya'], ya[index3])
        stim3['yb'] = np.append(stim3['yb'], yb[index3])
        stim3['xamc'] = np.append(stim3['xamc'], xamc[index3])
        stim3['yamc'] = np.append(stim3['yamc'], yamc[index3])
        stim3['yap'] = np.append(stim3['yap'], rtaph_yap(ya[index3],
                                                         yb[index3],
                                                         yamc[index3]))
        stim4['t'] = np.append(stim4['t'], t[index4])
        stim4['x'] = np.append(stim4['x'], x[index4])
        stim4['y'] = np.append(stim4['y'], y[index4])
        stim4['q'] = np.append(stim4['q'], q[index4])
        stim4['xa'] = np.append(stim4['xa'], xa[index4])
        stim4['xb'] = np.append(stim4['xb'], ya[index4])
        stim4['ya'] = np.append(stim4['ya'], ya[index4])
        stim4['yb'] = np.append(stim4['yb'], yb[index4])
        stim4['xamc'] = np.append(stim4['xamc'], xamc[index4])
        stim4['yamc'] = np.append(stim4['yamc'], yamc[index4])
        stim4['yap'] = np.append(stim4['yap'], rtaph_yap(ya[index4],
                                                         yb[index4],
                                                         yamc[index4]))

    print_inline("  Done.")

    return stim1, stim2, stim3, stim4
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def compute_stimstats(raw6file, band, eclipse):
    """
    Computes statistics for stim events for the post-CSP stim-based correction.

    :param raw6file: The name of the raw6 FITS file to read.

    :type raw6file: str

    :param band: The band to return the stim data for, either 'FUV' or 'NUV'.

    :type band: str

    :param eclipse: The eclipse number to return the stim data for.

    :type eclipse: int

    :returns: tuple -- Six-element tuple containing information on the trend
        in relative positions of the stims over time used for the Post-CSP
        stim correction.
    """

    print("Computing stim statistics and post-CSP corrections...")

    # Plate scale (in arcsec/mm).
    pltscl = 68.754932
    # Plate scale in arcsec per micron.
    aspum = pltscl/1000.0
    # [Future] If these are made module constants, can remove from this method.

    # Read in stim data from the FITS file.
    stim1, stim2, stim3, stim4 = raw6_to_stims(raw6file, band, eclipse)

    # Compute the mean positions (in arcseconds).
    stim1avg = [stim1['x'].mean()*aspum, stim1['y'].mean()*aspum]
    stim2avg = [stim2['x'].mean()*aspum, stim2['y'].mean()*aspum]
    stim3avg = [stim3['x'].mean()*aspum, stim3['y'].mean()*aspum]
    stim4avg = [stim4['x'].mean()*aspum, stim4['y'].mean()*aspum]

    print("Init: Number of stim photons:", len(stim1['t']),
          len(stim2['t']), len(stim3['t']), len(stim4['t']))
    print("Init: Mean x values at stim positions (arcsec):", stim1avg[0],
          stim2avg[0], stim3avg[0], stim4avg[0])
    print("Init: Mean x values at stim positions (arcsec):", stim1avg[1],
          stim2avg[1], stim3avg[1], stim4avg[1])
    print("Init: Mean y values at stim positions (micron):",
          stim1avg[1]/aspum, stim2avg[1]/aspum, stim3avg[1]/aspum,
          stim4avg[1]/aspum)

    # Compute the RMS around the mean (in arcseconds).
    stim1rms = [rms(stim1['x']*aspum), rms(stim1['y']*aspum)]
    stim2rms = [rms(stim2['x']*aspum), rms(stim2['y']*aspum)]
    stim3rms = [rms(stim3['x']*aspum), rms(stim3['y']*aspum)]
    stim4rms = [rms(stim4['x']*aspum), rms(stim4['y']*aspum)]

    # Compute the stim separation.
    stimsep = (((stim2avg[0]-stim1avg[0])+(stim4avg[0]-stim3avg[0])+
                (stim1avg[1]-stim3avg[1])+(stim2avg[1]-stim4avg[1]))/4.)
    print("Init: RMS  x values at stim positions (arcsec):", stim1rms[0],
          stim2rms[0], stim3rms[0], stim4rms[0])
    print("Init: RMS  y values at stim positions (arcsec):", stim1rms[1],
          stim2rms[1], stim3rms[1], stim4rms[1])
    print("Init: (arcsec): Stim sep =", stimsep, "    Average: X RMS =",
          (stim1rms[0]+stim2rms[0]+stim3rms[0]+stim4rms[0])/4.,
          "        Y RMS =",
          (stim1rms[1]+stim2rms[1]+stim3rms[1]+stim4rms[1])/4.)
    print("Raw stim separation is", stimsep)

    # Compute means and RMS values for each stim for each YA value stim1.
    for ya in range(32):
        ix = (stim1['ya'] == ya).nonzero()[0]
        ix = (stim2['ya'] == ya).nonzero()[0]
        ix = (stim3['ya'] == ya).nonzero()[0]
        ix = (stim4['ya'] == ya).nonzero()[0]

    # This returns the pre-CSP stim positions (because eclipse==0).
    avgstim = avg_stimpos(band, 0)

    # Compute Y scale and shift factors: yprime_as = (m * y_as) + B.
    y1, y2 = (stim1avg[1]+stim2avg[1])/2., (stim3avg[1]+stim4avg[1])/2.
    Y1, Y2 = ((avgstim['y1']+avgstim['y2'])/2.,
              (avgstim['y3']+avgstim['y4'])/2.)
    My = (Y1-Y2)/(y1-y2)
    By = (Y1-My*y1)/aspum
    print("Init: FODC: Y scale and shift (microns): My=", My, "By=", By)

    # Compute Y scale and shift factors: yprime_as = (m * y_as) + B.
    x1, x2 = (stim1avg[0]+stim3avg[0])/2., (stim2avg[0]+stim4avg[0])/2.
    X1, X2 = ((avgstim['x1']+avgstim['x3'])/2.,
              (avgstim['x2']+avgstim['x4'])/2.)
    Mx = (X1-X2)/(x1-x2)
    Bx = (X1-Mx*x1)/aspum
    print("Init: FODC: X scale and shift (microns): Mx=", Mx, "Bx=", Bx)

    stim1['xs'] = stim1['x']*Mx+Bx
    stim1['ys'] = stim1['y']*My+By
    stim2['xs'] = stim2['x']*Mx+Bx
    stim2['ys'] = stim2['y']*My+By
    stim3['xs'] = stim3['x']*Mx+Bx
    stim3['ys'] = stim3['y']*My+By
    stim4['xs'] = stim4['x']*Mx+Bx
    stim4['ys'] = stim4['y']*My+By

    # Compute the new mean positions (in arcseconds).
    stim1avgs = [stim1['xs'].mean()*aspum, stim1['ys'].mean()*aspum]
    stim2avgs = [stim2['xs'].mean()*aspum, stim2['ys'].mean()*aspum]
    stim3avgs = [stim3['xs'].mean()*aspum, stim3['ys'].mean()*aspum]
    stim4avgs = [stim4['xs'].mean()*aspum, stim4['ys'].mean()*aspum]

    print("Scal: Number of stim photons:", len(stim1['xs']),
          len(stim2['xs']), len(stim3['xs']), len(stim4['xs']))
    print("Scal: Mean x values at stim positions (arcsec):",
          stim1avgs[0], stim2avgs[0], stim3avgs[0], stim4avgs[0])
    print("Scal: Mean y values at stim positions (arcsec):",
          stim1avgs[1], stim2avgs[1], stim3avgs[1], stim4avgs[1])
    print("Scal: Mean y values at stim positions (microns):",
          stim1avgs[1]/aspum, stim2avgs[1]/aspum, stim3avgs[1]/aspum,
          stim4avgs[1]/aspum)

    # Compute the new RMS around the mean (in arcseconds).
    stim1rmss = [rms(stim1['xs']*aspum), rms(stim1['ys']*aspum)]
    stim2rmss = [rms(stim2['xs']*aspum), rms(stim2['ys']*aspum)]
    stim3rmss = [rms(stim3['xs']*aspum), rms(stim3['ys']*aspum)]
    stim4rmss = [rms(stim4['xs']*aspum), rms(stim4['ys']*aspum)]

    # Compute the stim separation.
    stimseps = (((stim2avgs[0]-stim1avgs[0])+(stim4avgs[0]-stim3avgs[0])+
                 (stim1avgs[1]-stim3avgs[1])+
                 (stim2avgs[1]-stim4avgs[1]))/4.)
    print("Scal: RMS  x values at stim positions (arcsec):",
          stim1rmss[0], stim2rmss[0], stim3rmss[0], stim4rmss[0])
    print("Init: RMS  y values at stim positions (arcsec):",
          stim1rmss[1], stim2rmss[1], stim3rmss[1], stim4rmss[1])
    print("Init: (arcsec): Stim sep =", stimseps, "   Average: X RMS =",
          (stim1rmss[0]+stim2rmss[0]+stim3rmss[0]+stim4rmss[0])/4.,
          "    Y RMS =",
          (stim1rmss[1]+stim2rmss[1]+stim3rmss[1]+stim4rmss[1])/4.)

    # Fit straight line to YA>2 and YB==2 points.
    # This could be written more efficiently, but it's an attempt at a faithful
    #  port of the GALEX code (from Tom Barlow) which was written in C.
    ix1 = ((stim1['ya'] > 2)&(stim1['yb'] == 2)).nonzero()[0]
    ix2 = ((stim2['ya'] > 2)&(stim2['yb'] == 2)).nonzero()[0]
    ix3 = ((stim3['ya'] > 2)&(stim3['yb'] == 2)).nonzero()[0]
    ix4 = ((stim4['ya'] > 2)&(stim4['yb'] == 2)).nonzero()[0]
    w8 = np.ones(len(ix1)+len(ix2)+len(ix3)+len(ix4))
    x8 = np.concatenate((stim1['yap'][ix1],
                         stim2['yap'][ix2],
                         stim3['yap'][ix3],
                         stim4['yap'][ix4]), axis=0)
    y8 = np.concatenate(
        (stim1['ys'][ix1]-stim1avgs[1]/aspum,
         stim2['ys'][ix2]-stim2avgs[1]/aspum,
         stim3['ys'][ix3]-stim3avgs[1]/aspum,
         stim4['ys'][ix4]-stim4avgs[1]/aspum), axis=0)
    print("NOTE: Found,", len(w8), "points for YA correction fit.")

    yac_coef1, yac_coef0 = np.polyfit(x8, y8, 1)

    print("Scal: YA correction coef for YB=2:", yac_coef0, yac_coef1)

    # Compute yb shift factors == zero for all.
    yac_ybs = np.zeros(8)
    coef0_yb = np.zeros(8)+yac_coef0
    coef1_yb = np.zeros(8)+yac_coef1

    # Set user slope adjustment. Use best slope adjustments from September 2010.
    # YB == 2...
    slope_scale = 1.04
    print("NOTE: Using slope scale of,", slope_scale, "for YB==2.")
    rr1 = yac_coef1*slope_scale
    rr0 = (yac_coef0 + (16.*yac_coef1))-(16.*rr1)
    coef0_yb[2] = rr0
    coef1_yb[2] = rr1
    print("New: YA correction coef (YB==2):", coef0_yb[2], coef1_yb[2])

    # YB == 3,4...
    slope_scale = 1.06
    print("NOTE: Using slope scale of,", slope_scale, "for YB==3.")
    rr1 = yac_coef1*slope_scale
    rr0 = (yac_coef0 + (16.*yac_coef1))-(16.*rr1)
    coef0_yb[3] = rr0
    coef1_yb[3] = rr1
    coef0_yb[4] = rr0
    coef1_yb[4] = rr1
    print("New: YA correction coef (YB==3):", coef0_yb[3], coef1_yb[3])
    print("NOTE: Using slope scale of,", slope_scale, "for YB==4.")
    print("New: YA correction coef (YB==4):", coef0_yb[4], coef1_yb[4])

    # Fill in look up array.
    yac = np.zeros([40, 8])
    for yb in range(8):
        for ya in range(40):
            yac[ya][yb] = ((coef0_yb[yb] +
                            (float(ya)*coef1_yb[yb])) +
                           yac_ybs[yb])

    stim1['yac'] = (yac[np.array(stim1['yap'], dtype='int64'),
                        np.array(stim1['yb'], dtype='int64')])
    stim2['yac'] = (yac[np.array(stim2['yap'], dtype='int64'),
                        np.array(stim2['yb'], dtype='int64')])
    stim3['yac'] = (yac[np.array(stim3['yap'], dtype='int64'),
                        np.array(stim3['yb'], dtype='int64')])
    stim4['yac'] = (yac[np.array(stim4['yap'], dtype='int64'),
                        np.array(stim4['yb'], dtype='int64')])

    # [Future] The section below could be re-written more elegantly.
    # [Future] Does this return the correct values for YB==1?
    for yb in range(8):
        ix = ((stim1['yb'] == yb)&(stim1['ya'] > 4)).nonzero()[0]
        s1m = ((stim1['ys']-stim1['yac'])[ix]*aspum).mean()
        s1r = rms((stim1['ys']-stim1['yac'])[ix]*aspum)
        if len(ix) > 0:
            print("Corrected stim 1: YB=", yb, " Num=", len(ix),
                  " Mean=", s1m, " RMS=", s1r)
    for yb in range(8):
        ix = ((stim2['yb'] == yb)&(stim2['ya'] > 4)).nonzero()[0]
        s2m = ((stim2['ys']-stim2['yac'])[ix]*aspum).mean()
        s2r = rms((stim2['ys']-stim2['yac'])[ix]*aspum)
        if len(ix) > 0:
            print("Corrected stim 2: YB=", yb, " Num=", len(ix),
                  " Mean=", s2m, " RMS=", s2r)
    for yb in range(8):
        ix = ((stim3['yb'] == yb)&(stim3['ya'] > 4)).nonzero()[0]
        s3m = ((stim3['ys']-stim3['yac'])[ix]*aspum).mean()
        s3r = rms((stim3['ys']-stim3['yac'])[ix]*aspum)
        if len(ix) > 0:
            print("Corrected stim 3: YB=", yb, " Num=", len(ix),
                  " Mean=", s3m, " RMS=", s3r)
    for yb in range(8):
        ix = ((stim4['yb'] == yb)&(stim4['ya'] > 4)).nonzero()[0]
        s4m = ((stim4['ys']-stim4['yac'])[ix]*aspum).mean()
        s4r = rms((stim4['ys']-stim4['yac'])[ix]*aspum)
        if len(ix) > 0:
            print("Corrected stim 4: YB=", yb, " Num=", len(ix),
                  " Mean=", s4m, " RMS=", s4r)

    return Mx, Bx, My, By, stimsep, yac
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def create_ssd(raw6file, band, eclipse, ssdfile=None):
    """
    Creates a Stim Separation Data (SSD) table file.

    :param raw6file: The name of the raw6 FITS file to read.

    :type raw6file: str

    :param band: The band to create the SSD file for, either 'FUV' or 'NUV'.

    :type band: str

    :param eclipse: The eclipse number to create the SSD file for.

    :type eclipse: int

    :param ssdfile: Name of stim separation data (SSD) output file to create.

    :type ssdfile: str

    :returns: tuple -- 2xN tuple containing slope and intercept of
        stim positions over time.
    """
    if ssdfile:
        print("Preparing SSD output file "+ssdfile)
        tbl = csv.writer(open(ssdfile, 'wb'), delimiter=' ',
                         quotechar='|', quoting=csv.QUOTE_MINIMAL)
        tbl.writerow(['|sct', '|stim_sep', '|stim_num', '|sep_fit'])

    # Plate scale (in arcsec/mm).
    pltscl = 68.754932
    # Plate scale in arcsec per micron.
    aspum = pltscl/1000.0
    # [Future] If these are made module constants, can remove from this method.

    stim1, stim2, stim3, stim4 = raw6_to_stims(raw6file, band, eclipse,
                                                                    margin=20.)
    # if not (stim1['t'].shape[0] and stim2['t'].shape[0] and
    #         stim3['t'].shape[0] and stim4['t'].shape[0]):
    #     # Missing a stim. Try again with the default margin of 90.
    #     stim1, stim2, stim3, stim4 = raw6_to_stims(raw6file, band, eclipse)
    #     if not (stim1['t'].shape[0] and stim2['t'].shape[0] and
    #             stim3['t'].shape[0] and stim4['t'].shape[0]):
    #         raise ValueError('Unable to locate a stim. Cant create SSD.')

    stimt = np.concatenate([stim1['t'], stim2['t'], stim3['t'], stim4['t']],
                           axis=0)
    sortt = np.argsort(stimt)
    stimt = stimt[sortt]
    stimix = np.concatenate([np.zeros(len(stim1['t']))+1,
                             np.zeros(len(stim2['t']))+2,
                             np.zeros(len(stim3['t']))+3,
                             np.zeros(len(stim4['t']))+4], axis=0)[sortt]
    stimx_as = (np.concatenate([stim1['x'], stim2['x'], stim3['x'],
                                stim4['x']], axis=0)*aspum)[sortt]
    stimy_as = (np.concatenate([stim1['y'], stim2['y'], stim3['y'],
                                stim4['y']], axis=0)*aspum)[sortt]
    pinc = 1000
    avt, sep, num = [], [], []

    for i in range(0, len(stimt)-pinc, pinc):
        ix1 = (stimix[i:i+pinc] == 1).nonzero()[0]
        ix2 = (stimix[i:i+pinc] == 2).nonzero()[0]
        ix3 = (stimix[i:i+pinc] == 3).nonzero()[0]
        ix4 = (stimix[i:i+pinc] == 4).nonzero()[0]
        sx1, sy1 = (
            np.mean(stimx_as[i:i+pinc][ix1]),
            np.mean(stimy_as[i:i+pinc][ix1])
            )
        sx2, sy2 = (
            np.mean(stimx_as[i:i+pinc][ix2]),
            np.mean(stimy_as[i:i+pinc][ix2])
            )
        sx3, sy3 = (
            np.mean(stimx_as[i:i+pinc][ix3]),
            np.mean(stimy_as[i:i+pinc][ix3])
            )
        sx4, sy4 = (
            np.mean(stimx_as[i:i+pinc][ix4]),
            np.mean(stimy_as[i:i+pinc][ix4])
            )
        stim_sep = ((sx2 - sx1) + (sx4 - sx3) + (sy1 - sy3) + (sy2 - sy4)) / 4.
        stim_avt = sum(stimt[i:i+pinc])/len(stimt[i:i+pinc])
        stim_num = len(ix1)+len(ix2)+len(ix3)+len(ix4)
        avt.append(stim_avt)
        sep.append(stim_sep)
        num.append(stim_num)

    m, C = np.polyfit(avt, sep, 1)
    fit = C+np.array(avt)*m

    if ssdfile:
        for i in range(len(avt)):
            tbl.writerow([avt[i], sep[i], num[i], fit[i]])

    return C, m
# ------------------------------------------------------------------------------
