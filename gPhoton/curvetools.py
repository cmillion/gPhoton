"""
.. module:: curvetools
   :synopsis: Functions for creation of lightcurves and components thereof.
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
from builtins import map
from builtins import str
from builtins import zip
import numpy as np
import pandas as pd
# gPhoton imports.
import gPhoton.cal as cal
import gPhoton.dbasetools as dbt
import gPhoton.galextools as gxt
import gPhoton.gQuery as gQuery
from gPhoton.gQuery import tscale
import gPhoton.MCUtils as mc
from gPhoton import __version__

# ------------------------------------------------------------------------------
def bg_contamination(band, skypos, radius, annulus=None):
    """
    Identify coadd MCAT sources that may contaminate the target observation.
    Uses a hard cutoff of 25 magnitude.

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: A two-element list containing the RA and DEC.

    :type skypos: list

    :param radius: Radius of the photometric aperture in degrees.

    :type radius: float

    :param annulus: A two-element list containing the inner and outer radius
        to use for background subtraction during aperture photometry, in
        degrees.

    :type annulus: list

    :returns: tuple -- Elements are: [1] number of coadd sources within the
        aperture, [2] number of coadd sources within the background annulus,
        [3] magnitude of brightnest source within the background annulus. If
        no annulus is specified then [2] and [3] are set to None.
    """

    bgsources = dbt.get_mags(band,skypos[0],skypos[1],
                    radius if annulus is None else annulus[1],25,verbose=0)
    if bgsources is None:
        return None, None, None
    mcatdist = mc.angularSeparation(skypos[0],skypos[1],
                                        bgsources['ra'],bgsources['dec'])
    apsourcecnt = np.where(mcatdist<=radius)[0].shape[0]
    bgsourcecnt,bgsourcemag=None,None
    if annulus is not None:
        bg_ix = np.where((annulus[0]<mcatdist) &
                               (mcatdist<=annulus[1]))
        bgsourcecnt = bg_ix[0].shape[0]
        if bgsourcecnt:
            bgsourcemag = bgsources[band]['mag'][bg_ix].max()
    return apsourcecnt, bgsourcecnt, bgsourcemag
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def gphot_params(band, skypos, radius, annulus=None, verbose=0, detsize=1.25,
                 stepsz=None, trange=None):
    """
    Populate a dict() with parameters that are constant over all bins.

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: A two-element list containing the RA and DEC.

    :type skypos: list

    :param radius: Radius of the photometric aperture in degrees.

    :type radius: float

    :param annulus: A two-element list containing the inner and outer radius
        to use for background subtraction during aperture photometry, in
        degrees.

    :type annulus: list

    :param verbose: Level of verbosity, 0 = minimum verbosity.

    :type verbose: int

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :param stepsz: Size of time bin to use, in seconds.

    :type stepsz: float

    :param trange: The start and end timestamps to consider, in GALEX time

    :type trange: list

    :returns: dict -- The set of parameters that are constant across all bins.
    """
    apsourcecnt, bgsourcecnt, bgsourcemag = bg_contamination(
                                        band, skypos, radius, annulus=annulus)

    return {'band':band, 'ra0':skypos[0], 'dec0':skypos[1], 'skypos':skypos,
            'trange':trange, 'radius':radius, 'annulus':annulus,
            'stepsz':stepsz, 'verbose':verbose, 'detsize':detsize,
            'apcorrect1':gxt.apcorrect1(radius, band),
            'apcorrect2':gxt.apcorrect2(radius, band),
            'detbg':gxt.detbg(mc.area(radius), band),
            'n_apersources':apsourcecnt,'n_bgsources':bgsourcecnt,
            'max_bgmag':bgsourcemag,'version':__version__,}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def xieta2colrow(xi, eta, band, detsize=1.25):
    """
    Convert detector xi and eta into col and row.

    :param xi: Sky-projected event "x" positions *in detetor coordinates*.

    :type xi: numpy.ndarray

    :param eta: Sky-projected event "y" positions *in detetor coordinates*.

    :type eta: numpy.ndarray

    :param band: The band that is being used, either 'FUV' or 'NUV'.

    :type band: str

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :returns: tuple -- A two-element tuple containing the column and row values.
    """

    flat, flatinfo = cal.flat(band)

    # Should be able to get npix from the header...
    npixx = flat.shape[0]
    npixy = flat.shape[1]

    pixsz = flatinfo['CDELT2']
    flatfill = detsize/(npixx*pixsz)

    col = (((xi/36000.)/(detsize/2.)*flatfill + 1.)/2. * npixx)
    row = (((eta/36000.)/(detsize/2.)*flatfill + 1.)/2. * npixy)

    return col, row
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def hashresponse(band, events, verbose=0):
    """
    Given detector xi, eta, return the response at each position.

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param events: The set of photon events.

    :type events: dict

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :returns: dict -- The photon event properties, updated with the response at
        each position.
    """

    # Hash out the response correction.
    if verbose:
        mc.print_inline("Applying the response correction.")
    flat, _ = cal.flat(band)

    events['col'], events['row'] = xieta2colrow(
        events['xi'], events['eta'], band)
    events['flat'] = flat[np.array(events['col'], dtype='int16'),
                          np.array(events['row'], dtype='int16')]
    events['scale'] = gxt.compute_flat_scale(events['t'], band)

    # [Future]: Separately do the binlinearly interpolated response.
    events['response'] = (events['flat']*events['scale'])

    return events
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def read_photons(photonfile, ra0, dec0, tranges, radius, verbose=0,
                 colnames=['t', 'x', 'y', 'xa', 'ya', 'q', 'xi', 'eta', 'ra',
                           'dec', 'flags']):
    """
    Read a photon list file and return a python dict() with the expected
        format.

    :param photonfile: Name of the photon event file to use.

    :type photonfile: str

    :param ra0: Right ascension of the targeted sky position, in degrees.

    :type ra0: float

    :param dec0: Declination of the targeted sky position, in degrees.

    :type dec0: float

    :param tranges: Set of time ranges from which to retrieve photon events,
        in GALEX time units

    :type tranges: list

    :param radius: The radius, in degrees, defining a cone on the sky that
        is centered on ra0 and dec0, from which to extract photons.

    :type radius: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param colnames: Labels of the columns found in the photon event file.

    :type colnames: list

    :returns: dict -- The set of photon events and their properties.
    """

    # [Future]: Consider moving this method to 'dbasetools'.

    if verbose:
        mc.print_inline('Reading photon list file: {f}'.format(f=photonfile))

    data = pd.io.parsers.read_csv(photonfile, names=colnames)
    ra, dec = np.array(data['ra']), np.array(data['dec'])
    angsep = mc.angularSeparation(ra0, dec0, ra, dec)

    ix = np.array([])
    for trange in tranges:
        cut = np.where((angsep <= radius) & (np.isfinite(angsep)))[0]
        ix = np.concatenate((ix, cut), axis=0)
    events = {'t':np.array(data['t'][ix])/tscale,
              'ra':np.array(data['ra'][ix]),
              'dec':np.array(data['dec'][ix]),
              'xi':np.array(data['xi'][ix]),
              'eta':np.array(data['eta'][ix]),
              'x':np.array(data['x'][ix]),
              'y':np.array(data['y'][ix])}

    return events
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def query_photons(band, ra0, dec0, tranges, radius, verbose=0, flag=0,
                  detsize=1.25):
    """
    Retrieve photons within an aperture from the database.

    :param band: Name of the band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: Right ascension of the targeted sky position, in degrees.

    :type ra0: float

    :param dec0: Declination of the targeted sky position, in degrees.

    :type dec0: float

    :param tranges: Set of time ranges from which to retrieve photon events,
        in GALEX time units

    :type tranges: list

    :param radius: The radius, in degrees, defining a cone on the sky that
        is centered on ra0 and dec0, from which to extract photons.

    :type radius: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param flag: Photon list flag value upon which to select. Default of 0
        corresponds to nominally corrected data with no issues. NOTE: 'Flag' is
        not a reliable way to parse data at this time. You should compare
        timestamps against the aspect file.

    :type flag: int

    :returns: dict -- The set of photon events with their properties.
    """

    # [Future]: This should probably be moved to 'dbasetools'.
    stream = []
    if verbose:
        mc.print_inline(
            "Retrieving photons within {rad} degrees of [{r}, {d}]".format(
                rad=radius, r=ra0, d=dec0))
    for trange in tranges:
        if verbose:
            mc.print_inline(" and between {t0} and {t1}.".format(t0=trange[0],
                                                                 t1=trange[1]))

        # [Future]: This call to fGetTimeRanges prevents the downloading of
        # events with bad aspect solutions which currently have incorrect
        # quality flags (of zero) in the photon database.
        trs = dbt.fGetTimeRanges(band, [ra0, dec0], trange=trange,
                                 detsize=detsize)
        if not trs.any():
            continue
        for tr in trs:
            thisstream = gQuery.getArray(
                gQuery.allphotons(band, ra0, dec0, tr[0], tr[1], radius,
                                  flag=flag), verbose=verbose, retries=100)
            stream.extend(thisstream)

    stream = np.array(stream, 'f8').T
    colnames = ['t', 'ra', 'dec', 'xi', 'eta', 'x', 'y', 'flag', 'q']
    dtypes = ['f8', 'f8', 'f8', 'f4', 'f4', 'f4', 'f4', 'i4', 'f4']
    try:
        cols = list(map(np.asarray, stream, dtypes))
        events = dict(list(zip(colnames, cols)))
        # Adjust the timestamp by tscale.
        events['t'] /= tscale
        return events
    except TypeError:
        return None
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def pullphotons(band, ra0, dec0, tranges, radius, verbose=0, flag=0,
                photonfile=None, detsize=1.25):
    """
    Reads photon list data from the MAST database using a cone search.

    :param band: Name of the band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: Right ascension of the targeted sky position, in degrees.

    :type ra0: float

    :param dec0: Declination of the targeted sky position, in degrees.

    :type dec0: float

    :param tranges: Set of time ranges from which to retrieve photon events,
        in GALEX time units

    :type tranges: list

    :param radius: The radius, in degrees, defining a cone on the sky that
        is centered on ra0 and dec0, from which to extract photons.

    :type radius: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param photonfile: Name of photon event file to use (if reading from disk).

    :type photonfile: str

    :param flag: Photon list flag value upon which to select. Default of 0
        corresponds to nominally corrected data with no issues. NOTE: 'Flag' is
        not a reliable way to parse data at this time. You should compare
        timestamps against the aspect file.

    :type flag: int

    :returns: dict -- The set of photon events and their properties.
    """

    if photonfile:
        events = read_photons(photonfile, ra0, dec0, tranges, radius,
                              verbose=verbose)
    else:
        events = query_photons(band, ra0, dec0, tranges, radius,
                               verbose=verbose, flag=flag, detsize=detsize)

    try:
        events = hashresponse(band, events, verbose=verbose)
        return events
    except TypeError:
        return None
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def aperture_error(counts, expt, bgcounts=0):
    """
    The estimated error in the countrate within the aperture, by adding
        together the counting error within the aperture and background (if
        provided) in quadrature.

    :param counts: Total counts within the aperture.

    :type counts: int

    :param expt: The exposure time in seconds.

    :type extp: float

    :param bgcounts: The total background counts within the aperture.

    :type bgcounts: int

    :returns: float -- The error in the counts within the aperture.
    """

    return np.sqrt(counts/expt**2+bgcounts/expt**2)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def reduce_lcurve(bin_ix, region_ix, data, function, dtype='float64'):
    """
    Produces light curve columns by iteratively applying 'function' to 'data'
        within 'region_ix' over 'bin_ix'.

    :param bin_ix: Array indices designating which events are in the time bin
        of interest.

    :type bin_ix: numpy.ndarray

    :param region_ix: Array indices designating which events are in the spatial
        region of interest (e.g. the photometric aperture).

    :type region_ix: numpy.ndarray

    :param data: The data to apply the function on.

    :type data: numpy.ndarray

    :param function: The function to apply to the data.

    :type function: function

    :param dtype: The data type.

    :type dtype: str

    :returns: numpy.ndarray -- The light curve columns.
    """

    bin_num = np.unique(bin_ix)
    output = np.empty(len(bin_num))

    for i, b in enumerate(bin_num):
        if len(np.where(bin_ix[region_ix]==b)[0])==0:
            continue
        try:
            ix = region_ix[0][np.where(bin_ix[region_ix] == b)]
            output[i] = function(data[ix])
        except ValueError:
            output[i] = np.nan
        except IndexError:
            output[i] = np.nan
        except:
            raise

    return np.array(output, dtype=dtype)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def recoverywarning(band, bin_ix, events, verbose=0):
    """
    Test whether the bin contains data that was collected during a spacecraft
        recovery period (e.g. FUV cycling) as defined by the lookup table in
        galextools.recovery_tranges().

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param bin_ix: Array indices designating which events are in the time bin
        of interest.

    :type bin_ix: numpy.ndarray

    :param events: Set of photon events to check if they are near a masked
        detector region.

    :type events: dict

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int
    """
    tranges = gxt.recovery_tranges()
    for trange in tranges:
        t = np.array(events['photons']['t'])[bin_ix]
        if ((trange[0] <= t) & (trange[1] >= t)).any():
            return True
    return False
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def caiwarning(band, bin_ix, events, verbose=0):
    """
    Test whether a bin contains data from the first 3 legs of an FUV
        observation as part of the calibration (CAI) survey.

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param bin_ix: Array indices designating which events are in the time bin
        of interest.

    :type bin_ix: numpy.ndarray

    :param events: Set of photon events to check if they are near a masked
        detector region.

    :type events: dict

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int
    """
    if band == 'FUV':
        for trange in dbt.distinct_tranges(np.array(events['photons']['t'])[bin_ix]):
            t = np.median(trange)
            obsdata = gQuery.getArray(gQuery.obstype_from_t(t))
            if not obsdata:
                continue
            if (str(obsdata[0][0]) is 'CAI') and (obsdata[0][5] <= 3):
                return True
    return False
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def maskwarning(band, bin_ix, events, verbose=0, mapkey='H', mode=None):
    """
    Test if any given events are near a masked detector region.

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param bin_ix: Array indices designating which events are in the time bin
        of interest.

    :type bin_ix: numpy.ndarray

    :param events: Set of photon events to check if they are near a masked
        detector region.

    :type events: dict

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param mapkey: Text code indicating whether to use the hotspot mask ("H")
        or the flat (edge) mask ("E").

    :type mapkey: str

    :returns: bool -- Returns True/False whether a given set of events are near
        a masked detector region.
    """

    maps = {'H':cal.mask, 'E':cal.flat}

    img, _ = maps[mapkey](band, buffer=True)

    if mode is None:
        reg_ix = np.where(events['photons']['col'][bin_ix]) # i.e. all of them
    elif mode is 'aper':
        reg_ix = np.where(
            mc.angularSeparation(
                events['params']['skypos'][0],
                events['params']['skypos'][1],
                events['photons']['ra'],
                events['photons']['dec'])[bin_ix] <= events['params']['radius'])
    elif mode is 'bg':
        if not events['params']['annulus']:
            return False
        reg_ix = np.where(
            (mc.angularSeparation(
                events['params']['skypos'][0],
                events['params']['skypos'][1],
                events['photons']['ra'],
                events['photons']['dec'])[bin_ix] <= (
                    events['params']['annulus'][0])) &
            (mc.angularSeparation(
                events['params']['skypos'][0],
                events['params']['skypos'][1],
                events['photons']['ra'],
                events['photons']['dec'])[bin_ix] < (
                    events['params']['annulus'][1])))
    else:
        print('Unknown mask flag mode of: {m}'.format(m=mode))
        raise ValueError("Unknown mask flag mode.")

    for xoff in [-1, 0, 1]:
        for yoff in [-1, 0, 1]:
            if np.shape(np.where(
                    img[np.array(
                        events['photons']['col'][bin_ix][reg_ix],
                        dtype='int32')+xoff,
                        np.array(
                            events['photons']['row'][bin_ix][reg_ix],
                            dtype='int32')+yoff] == 0))[1] > 0:
                return True

    return False#True if len(ix[0]) else False
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def lowresponsewarning(bin_ix, events, verbose=0, ratio=0.7):
    """
    Checks for anomalously low response values in the data of interest, which
        could indicate data on poorly characterized or behaved detector regions.

    :param bin_ix: Array indices designating which events are in the time bin
        of interest.

    :type bin_ix: numpy.ndarray

    :param events: Set of photon events to check if there is a low response.

    :type events: dict

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param ratio: The value that defines a low response, between 0 and 1.
        (Where a response value of 1 indicates a "perfect" response value w/ no
        correction.)

    :type ratio: float

    :returns: bool -- Returns True/False whether a given set of events contain
        any on a low response region of the detector.
    """

    ix = np.where(events['photons']['response'][bin_ix] < 0.7)

    return True if len(ix[0]) else False
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def exptimewarning(bin_ix, events, verbose=0, ratio=0.5):
    """
    Passes a warning if the effective exposure time within a bin is
        significantly less than the raw exposure time, which might produce
        anomalous values due to counting statistics or be a symptom of a problem
        in the exposure time correction for this bin.

    :param bin_ix: Array indices designating which events are in the time bin
        of interest.

    :type bin_ix: numpy.ndarray

    :param events: Set of photon events to check if there is an effective
        exposure time warning.

    :type events: dict

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param ratio: The ratio of effective to raw exposure time in a bin below
        which the bin will be flagged.

    :type ratio: float

    :returns: bool -- Returns True/False whether a bin has a low exposure.
    """

    return (events['exptime'][bin_ix]/
            (events['t1'][bin_ix]-events['t0'][bin_ix]) < ratio)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def nonlinearitywarning(band, bin_ix, events, verbose=0):
    """
    Flag count rates above the 10% local nonlinearty dropoff, per the
        calibration paper.

    :param band: The band that is being used, either 'FUV' or 'NUV'.

    :type band: str

    :param bin_ix: Array indices designating which events are in the time bin
        of interest.

    :type bin_ix: numpy.ndarray

    :param events: Set of photon events to check if they are in the
        non-linearity regime.

    :type events: dict

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :returns: bool -- Returns True/False whether a given set of events are at
        the non-linearity regime.
    """

    cps_10p_rolloff = {'NUV':311, 'FUV':109}

    cps = events['flat_counts'][bin_ix]/events['exptime'][bin_ix]

    return True if cps >= cps_10p_rolloff[band] else False
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def detedgewarning(bin_ix, events, verbose=0, valid_detrad=0.5):
    """
    Assigns warning flags if any of the events of interest are adjacent
        to the detector edge as defined by a radius of valid_detrad in degrees.

    :param bin_ix: Array indices designating which events are in the time bin
        of interest.

    :type bin_ix: numpy.ndarray

    :param events: Set of photon events to check if they are near the detector
        edge.

    :type events: dict

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param valid_detrad: The radius, in degrees, beyond which an edge warning is
        raised.

    :type valid_detrad: float

    :returns: bool -- Returns True/False whether a given set of events are too
        close to the edge of the detector.
    """

    ix = np.where(mc.distance(events['photons']['col'][bin_ix],
                              events['photons']['row'][bin_ix], 400, 400)*
                  gxt.aper2deg(4) >= valid_detrad)

    return True if len(ix[0]) else False
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def getflags(band, bin_ix, events, verbose=0):
    """
    Pass flags if data meets conditions that are likely to create
        misleading photometry. The flags are binary, with bins set as follows:
        1 - 'hotspot' - aperture events in pixels contiguous to a masked hotspot
        2 - 'mask edge' - aperture events in pixels contiguous to the detector
        edge
        4 - 'exptime' - bin contains < 50% exposure time coverage
        8 - 'respose' - events weighted with response < 0.7
        16 - 'nonlinearity' - local countrate exceeds 10% response dropoff
        32 - 'detector edge' - events outside of 0.5 degrees of detector center
        64 - 'bg hotspot' - annulus events in pixels contiguous to a masked
        hotspot
        128 - 'bg mask' - annulus events in pixels contiguous to detector edge

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param bin_ix: Array indices designating which events are in the time bin
        of interest.

    :type bin_ix: numpy.ndarray

    :param events: Set of photon events to check for warning flags.

    :type events: dict

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :returns: numpy.ndarray -- The array of flags for each photon event.
    """

    bin_num = np.unique(bin_ix)
    flags = np.zeros(len(bin_num))

    for i, b in enumerate(bin_num):
        ix = np.where(bin_ix == b)
        if len(ix):
            #ix = bin_ix[np.where(bin_ix == bin)]
            try:
                if maskwarning(band, ix, events, mapkey='H',
                               mode='aper', verbose=verbose):
                    flags[i] += 1
                if maskwarning(band, ix, events, mapkey='E',
                               mode='aper', verbose=verbose):
                    flags[i] += 2
                if exptimewarning(i, events, verbose=verbose):
                    flags[i] += 4
                if lowresponsewarning(ix, events, verbose=verbose):
                    flags[i] += 8
                if nonlinearitywarning(band, i, events, verbose=verbose):
                    flags[i] += 16
                if detedgewarning(ix, events, verbose=verbose):
                    flags[i] += 32
                if maskwarning(band, ix, events, mapkey='H',
                               mode='bg', verbose=verbose):
                    flags[i] += 64
                if maskwarning(band, ix, events, mapkey='E',
                               mode='bg', verbose=verbose):
                    flags[i] += 128
                #if caiwarning(band, ix, events, verbose=verbose):
                #    flags[i] += 256
                if recoverywarning(band, ix, events, verbose=verbose):
                    flags[i] += 512
            except:
                raise
        else:
            flags[i] = np.nan

    return np.array(flags)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def quickmag(band, ra0, dec0, tranges, radius, annulus=None, stepsz=None,
             verbose=0, detsize=1.25, coadd=False):
    """
    Primary wrapper function for generating and synthesizing all of the
        parameters and calculations necessary to create light curves.

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: Right ascension, in degrees, of the target position.

    :type ra0: float

    :param dec0: Declination, in degrees, of the target position.

    :type dec0: float

    :param tranges: Set of time ranges to query within in GALEX time seconds.

    :type tranges: list

    :param radius: The radius of the  photometric aperture, in degrees.

    :type radius: float

    :param annulus: Radii of the inner and outer extents of the background
        annulus, in degrees.

    :type annulus: list

    :param stepsz: The size of the time bins to use, in seconds.

    :type stepsz: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :param coadd: Set to True if calculating a total flux instead of flux
        from each time bin.

    :type coadd: bool

    :returns: dict -- The light curve, including input parameters.
    """

    if verbose:
        mc.print_inline("Retrieving all of the target events.")
    searchradius = radius if annulus is None else annulus[1]
    data = pullphotons(band, ra0, dec0, tranges, searchradius, verbose=verbose,
                       detsize=detsize)
    if not data:
        return None

    if verbose:
        mc.print_inline("Binning data according to requested depth.")

    # Multiple ways of defining bins
    try:
        trange = [np.array(tranges).min(), np.array(tranges).max()]
    except ValueError:
        trange = tranges

    if coadd:
        bins = np.array(trange)
    elif stepsz:
        bins = np.append(np.arange(trange[0], trange[1], stepsz), max(trange))
    else:
        bins = np.unique(np.array(tranges).flatten())

    lcurve = {'params':gphot_params(band, [ra0, dec0], radius, annulus=annulus,
                                    verbose=verbose, detsize=detsize,
                                    stepsz=stepsz, trange=trange)}

    # This is equivalent in function to np.digitize(data['t'],bins) except
    # that it's much, much faster. See numpy issue #2656 at
    # https://github.com/numpy/numpy/issues/2656
    bin_ix = np.searchsorted(bins, data['t'], "right")
    try:
        lcurve['t0'] = bins[np.unique(bin_ix)-1]
        lcurve['t1'] = bins[np.unique(bin_ix)]
        lcurve['exptime'] = np.array(
            dbt.compute_exptime(band,
                                tranges if coadd else list(zip(
                                    lcurve['t0'], lcurve['t1'])),
                                verbose=verbose, coadd=coadd, detsize=detsize,
                                skypos=[ra0, dec0]))
    except IndexError:
        if np.isnan(data['t']):
            if verbose:
                mc.print_inline(
                    "No valid data available in {t}".format(t=tranges))
        lcurve['t0'] = np.array([np.nan])
        lcurve['t1'] = np.array([np.nan])
        lcurve['exptime'] = np.array([0])

    angSep = mc.angularSeparation(ra0, dec0, data['ra'], data['dec'])

    aper_ix = np.where(angSep <= radius)
    lcurve['t0_data'] = reduce_lcurve(bin_ix, aper_ix, data['t'], np.min)
    lcurve['t1_data'] = reduce_lcurve(bin_ix, aper_ix, data['t'], np.max)
    lcurve['t_mean'] = reduce_lcurve(bin_ix, aper_ix, data['t'], np.mean)
    lcurve['q_mean'] = reduce_lcurve(bin_ix, aper_ix, data['q'], np.mean)
    lcurve['counts'] = reduce_lcurve(bin_ix, aper_ix, data['t'], len)
    lcurve['flat_counts'] = reduce_lcurve(bin_ix, aper_ix,
                                          1./data['response'], np.sum)
    lcurve['responses'] = reduce_lcurve(bin_ix, aper_ix, data['response'],
                                        np.mean)
    lcurve['detxs'] = reduce_lcurve(bin_ix, aper_ix, data['col'], np.mean)
    lcurve['detys'] = reduce_lcurve(bin_ix, aper_ix, data['row'], np.mean)
    lcurve['detrad'] = mc.distance(lcurve['detxs'], lcurve['detys'], 400, 400)
    lcurve['racent'] = reduce_lcurve(bin_ix, aper_ix, data['ra'], np.mean)
    lcurve['deccent'] = reduce_lcurve(bin_ix, aper_ix, data['dec'], np.mean)

    skybgmcatdata = dbt.get_mcat_data([ra0, dec0], radius)
    lcurve['mcat_bg'] = lcurve['exptime']*np.array(
        [dbt.mcat_skybg(band, [ra0, dec0], radius, trange=tr,
                        mcat=skybgmcatdata, verbose=verbose)
         for tr in zip(lcurve['t0'], lcurve['t1'])])

    if annulus is not None:
        annu_ix = np.where((angSep > annulus[0]) & (angSep <= annulus[1]))
        lcurve['bg_counts'] = reduce_lcurve(bin_ix, annu_ix, data['t'], len)
        lcurve['bg_flat_counts'] = reduce_lcurve(
            bin_ix, annu_ix, data['response'], np.sum)
        lcurve['bg'] = (mc.area(radius)*lcurve['bg_flat_counts'] /
                        (mc.area(annulus[1])-mc.area(annulus[0])))
    else:
        lcurve['bg_counts'] = np.zeros(len(lcurve['counts']))
        lcurve['bg_flat_counts'] = np.zeros(len(lcurve['counts']))
        lcurve['bg'] = np.zeros(len(lcurve['counts']))

    lcurve['cps'] = lcurve['flat_counts']/lcurve['exptime']
    lcurve['cps_err'] = aperture_error(lcurve['flat_counts'], lcurve['exptime'])
    lcurve['cps_bgsub'] = (lcurve['flat_counts']-
                           lcurve['bg'])/lcurve['exptime']
    lcurve['cps_bgsub_err'] = aperture_error(
        lcurve['flat_counts'], lcurve['exptime'], bgcounts=lcurve['bg'])
    lcurve['cps_mcatbgsub'] = (lcurve['flat_counts']-
                               lcurve['mcat_bg'])/lcurve['exptime']
    lcurve['cps_mcatbgsub_err'] = aperture_error(
        lcurve['flat_counts'], lcurve['exptime'], bgcounts=lcurve['mcat_bg'])
    lcurve['flux'] = gxt.counts2flux(lcurve['cps'], band)
    lcurve['flux_err'] = gxt.counts2flux(lcurve['cps_err'], band)
    lcurve['flux_bgsub'] = gxt.counts2flux(lcurve['cps_bgsub'], band)
    lcurve['flux_bgsub_err'] = gxt.counts2flux(lcurve['cps_bgsub_err'], band)
    lcurve['flux_mcatbgsub'] = gxt.counts2flux(lcurve['cps_mcatbgsub'], band)
    lcurve['flux_mcatbgsub_err'] = gxt.counts2flux(
        lcurve['cps_mcatbgsub_err'], band)

    # NOTE: These conversions to mag can throw logarithm warnings if the
    # background is brighter than the source, resuling in a negative cps which
    # then gets propagated as a magnitude of NaN.
    lcurve['mag'] = gxt.counts2mag(lcurve['cps'], band)
    lcurve['mag_err_1'] = (lcurve['mag'] -
                           gxt.counts2mag(lcurve['cps'] + lcurve['cps_err'],
                                          band))
    lcurve['mag_err_2'] = (gxt.counts2mag(lcurve['cps'] -
                                          lcurve['cps_err'], band) -
                           lcurve['mag'])
    lcurve['mag_bgsub'] = gxt.counts2mag(lcurve['cps_bgsub'], band)
    lcurve['mag_bgsub_err_1'] = (lcurve['mag_bgsub'] -
                                 gxt.counts2mag(lcurve['cps_bgsub'] +
                                                lcurve['cps_bgsub_err'], band))
    lcurve['mag_bgsub_err_2'] = (gxt.counts2mag(lcurve['cps_bgsub'] -
                                                lcurve['cps_bgsub_err'], band) -
                                 lcurve['mag_bgsub'])
    lcurve['mag_mcatbgsub'] = gxt.counts2mag(lcurve['cps_mcatbgsub'], band)
    lcurve['mag_mcatbgsub_err_1'] = (lcurve['mag_mcatbgsub'] -
                                     gxt.counts2mag(lcurve['cps_mcatbgsub'] +
                                                    lcurve['cps_mcatbgsub_err'],
                                                    band))
    lcurve['mag_mcatbgsub_err_2'] = (gxt.counts2mag(lcurve['cps_mcatbgsub'] -
                                                    lcurve['cps_mcatbgsub_err'],
                                                    band) -
                                     lcurve['mag_mcatbgsub'])

    lcurve['photons'] = data

    lcurve['flags'] = getflags(band, bin_ix, lcurve, verbose=verbose)

    return lcurve
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def get_curve(band, ra0, dec0, radius, annulus=None, stepsz=None,
              trange=None, tranges=None, verbose=0, coadd=False, minexp=1.,
              maxgap=1., detsize=1.1):
    """
    Wraps quickmag() to make it ensure some proper parameter formatting and
        therefore make it slightly more user friendly.

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: Right ascension, in degrees, of the target position.

    :type ra0: float

    :param dec0: Declination, in degrees, of the target position.

    :type dec0: float

    :param radius: The radius of the photometric aperture, in degrees.

    :type radius: float

    :param annulus: Radii of the inner and outer extents of the background
        annulus, in degrees.

    :type annulus: list

    :param stepsz: The size (depth) of the time bins to use, in seconds.

    :type stepsz: float

    :param trange: Minimum and maximum time range to make a light curve,
        in GALEX time seconds.

    :type trange: list

    :param tranges: Set of time ranges to query within in GALEX time seconds.

    :type tranges: list

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param coadd: Set to True if calculating a total flux instead of flux
        from each time bin.

    :type coadd: bool

    :param minexp: Minimum gap size, in seconds, for data to be considered
        contiguous.

    :type minexp: float

    :param maxgap: Maximum gap size, in seconds, for data to be considered
        contiguous.

    :type maxgap: float

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :returns: dict -- The light curve, including input parameters.
    """

    # Temporary # HACK:
    if trange is not None:
        trange[0]-=trange[0]%0.005
        trange[1]-=trange[1]%0.005

    skyrange = [np.array(annulus).max().tolist() if annulus else radius,
                np.array(annulus).max().tolist() if annulus else radius,]
    if verbose:
        mc.print_inline("Getting exposure ranges.")
    if tranges is None:
        tranges = dbt.fGetTimeRanges(band, [ra0, dec0], trange=trange,
                                     maxgap=maxgap, minexp=minexp,
                                     verbose=verbose, detsize=detsize)
    if not len(np.array(tranges).flatten()):
        if verbose:
            mc.print_inline(
                "No exposure time at this location: [{ra},{dec}]".format(
                    ra=ra0, dec=dec0))
        return None
    if trange is not None and stepsz:
        tranges = [trange] # affix the time ranges under this condition
                           # for reproducibility and bin-matching
    lcurve = quickmag(band, ra0, dec0, tranges, radius, annulus=annulus,
                      stepsz=stepsz, verbose=verbose, coadd=coadd)
    if not lcurve:
        return None

    if verbose:
        mc.print_inline("Done.")
        mc.print_inline("")

    return lcurve
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def write_curve(band, ra0, dec0, radius, csvfile=None, annulus=None,
                stepsz=None, trange=None, tranges=None, verbose=0, coadd=False,
                iocode='w', detsize=1.1, overwrite=False, minexp=1., maxgap=1.,
                minimal_output=False, photoncsvfile=None, addhdr=False,
                commentchar='|'):
    """
    Generates a lightcurve and optionally writes the data to a CSV file.

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: Right ascension, in degrees, of the target position.

    :type ra0: float

    :param dec0: Declination, in degrees, of the target position.

    :type dec0: float

    :param radius: The radius of the photometric aperture, in degrees.

    :type radius: float

    :param csvfile: Name of the photon event CSV file to use for the lightcurve.

    :type csvfile: str

    :param annulus: Radii of the inner and outer extents of the background
        annulus, in degrees.

    :type annulus: list

    :param stepsz: The size of the time bins to use, in seconds.

    :type stepsz: float

    :param trange: Minimum and maximum timew within which to make a lightcurve,
        in GALEX time seconds.

    :type trange: list

    :param tranges: Set of time ranges to query within in GALEX time seconds.

    :type tranges: list

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param coadd: Set to True if calculating a total flux instead of flux
        from each time bin.

    :type coadd: bool

    :param iocode: The code to use when writing the output file.

    :type iocode: str

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :param overwrite: If True, overwite an existing output file.

    :type overwrite: bool

    :param minexp: Minimum gap size, in seconds, for data to be considered
        contiguous.

    :type minexp: float

    :param maxgap: Maximum gap size, in seconds, for data to be considered
        contiguous.

    :type maxgap: float

    :param minimal_output: If True, produce an output file with a minimum
        number of columns.

    :type minimal_output: bool

    :param photoncsvfile: Name of the photon event CSV file to write the
        photon list data to.

    :type photoncsvfile: str

    :returns: dict -- The light curve, including input parameters.
    """

    data = get_curve(band, ra0, dec0, radius, annulus=annulus, stepsz=stepsz,
                     trange=trange, tranges=tranges, verbose=verbose,
                     coadd=coadd, minexp=minexp, maxgap=maxgap,
                     detsize=detsize)
    if not data:
        if verbose:
            print('No events available at the requested location and time(s).')
        return None

    exclude_keys = ['photons', 'params']

    if minimal_output:
        # Initialize a list that will define the columns included in the
        # minimal output option.
        output_columns = ['t_mean']

        # Only include the annulus-corrected fluxes if an annulus was specified.
        if annulus:
            output_columns += ['flux_bgsub', 'flux_bgsub_err']

            # Then include these columns in all cases.
            output_columns += ['flux_mcatbgsub', 'flux_mcatbgsub_err',
                               'flags', 'counts', 'exptime', 'detrad']
    else:
        # Otherwise, define the order of the output columns to match the
        # User Guide sorting.
        time_related_cols = ['t0', 't1', 't_mean', 't0_data', 't1_data']
        annulus_bkg_corrected_cols = ['cps_bgsub', 'cps_bgsub_err',
                                      'flux_bgsub', 'flux_bgsub_err',
                                      'mag_bgsub', 'mag_bgsub_err_1',
                                      'mag_bgsub_err_2']
        mcat_bkg_corrected_cols = ['cps_mcatbgsub', 'cps_mcatbgsub_err',
                                   'flux_mcatbgsub', 'flux_mcatbgsub_err',
                                   'mag_mcatbgsub', 'mag_mcatbgsub_err_1',
                                   'mag_mcatbgsub_err_2']
        bkg_uncorrected_cols = ['cps', 'cps_err',
                                'flux', 'flux_err',
                                'mag', 'mag_err_1', 'mag_err_2']
        total_counts_cols = ['counts', 'flat_counts', 'bg_counts',
                             'bg_flat_counts']
        calibration_cols = ['exptime', 'bg', 'mcat_bg', 'responses', 'detxs',
                            'detys', 'detrad', 'racent', 'deccent', 'q_mean', 'flags']
        output_columns = (time_related_cols + annulus_bkg_corrected_cols +
                          mcat_bkg_corrected_cols + bkg_uncorrected_cols +
                          total_counts_cols + calibration_cols)

    # The output columns defined here must match all of those in the return
    # data dict, if not, raise an error so it can be fixed by developers.
    if minimal_output:
        # Make sure each of these are included in data.keys()
        if not set(output_columns) <= set(
                [x for x in list(data.keys()) if x not in exclude_keys]):
            raise ValueError("Some output columns are not present in the"
                             " data structure.  A developer needs to fix this.")
    else:
        if set(output_columns) != set(
                [x for x in list(data.keys()) if x not in exclude_keys]):
            raise ValueError("Output columns do not match those returned in"
                             " data structure.  A developer needs to fix this.")

    if csvfile:
        if verbose:
            mc.print_inline('Building output data frame.')

        frame, columns = {}, []

        for k in output_columns:
            if k in exclude_keys:
                continue
            frame[k] = data[k]
            columns += [k]

        try:
            output = pd.DataFrame(frame)
        except:
            if verbose > 1:
                print('Unable to build dataframe.')
            raise ValueError("Unable to build dataframe.")
        try:
            if addhdr:
                with open(csvfile, iocode) as f:
                    for k in data['params'].keys():
                        f.write('{c} {k} = {v}\n'.format(
                            c=commentchar,k=k,v=data['params'][k]))
            output.to_csv(csvfile, index=False, mode='a' if addhdr else iocode,
                          columns=columns)
        except:
            print('Unable to write to: '+str(csvfile))
    else:
        if verbose > 2:
            print("No CSV file requested.")

        if verbose or (not verbose and not csvfile):
            print("AB Magnitudes:               ")
            print(data['mag'])

    if photoncsvfile:
        pd.DataFrame(data['photons']).to_csv(photoncsvfile, index=False)

    return data
# ------------------------------------------------------------------------------
