"""
.. module:: dbasetools
   :synopsis: Contains tools for working with data from the database that are
       used by a number of different modules.

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
from builtins import range
from builtins import str
from builtins import zip
import numpy as np
# gPhoton imports.
from gPhoton.galextools import GPSSECS, zpmag, aper2deg
import gPhoton.gQuery as gQuery
from gPhoton.gQuery import tscale
from gPhoton.MCUtils import print_inline, area, angularSeparation

# ------------------------------------------------------------------------------
def get_aspect(band, skypos, trange=[6e8, 11e8], verbose=0, detsize=1.25):
    """
    Get aspect solution in a dict() for the given time range.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param trange: Minimum and maximum time (in GALEX time) to consider.

    :type trange: list

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param detsize: The effective detector diameter, in degrees.

    :type detsize: float

    :returns: dict -- The aspect solution parameters.
    """

    asp = np.array(gQuery.getArray(gQuery.aspect_skypos(skypos[0], skypos[1],
                                                        detsize=detsize),
                                   verbose=verbose))

    data = {'eclipse':np.array(asp[:, 0], dtype='int16'), 'filename':asp[:, 1],
            't':np.array(asp[:, 2], dtype='float64')/tscale,
            'ra':np.array(asp[:, 3], dtype='float64'),
            'dec':np.array(asp[:, 4], dtype='float64'),
            'twist':np.array(asp[:, 5], dtype='float64'),
            'flag':np.array(asp[:, 6], dtype='int8'),
            'ra0':np.array(asp[:, 7], dtype='float64'),
            'dec0':np.array(asp[:, 8], dtype='float64'),
            'twist0':np.array(asp[:, 9], dtype='float64')}

    ix = np.where((data['t'] > trange[0]) & (data['t'] < trange[1]) &
                  (angularSeparation(skypos[0], skypos[1],
                                     data['ra'], data['dec']) <= detsize/2.))

    for key in list(data.keys()):
        data[key] = data[key][ix]

    return data
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def distinct_tranges(times, maxgap=1.):
    """
    Produces a list of pairs of start / stop times delimiting distinct
        unique time ranges, given that gaps of >maxgap initiate a new time
        period.

    :param times: A set of time stamps to extract unique time ranges from.

    :type times: list

    :param maxgap:  Maximum gap size, in seconds, for data to be considered
        contiguous.

    :type maxgap: float

    :returns: list -- The set of distinct time ranges.
    """

    ix = np.where(times[1:]-times[:-1] > maxgap)

    ixs = [-1] + list(ix[0]) + [len(times)-1]

    return [[times[ixs[i]+1], times[ixs[i+1]]] for i in range(len(ixs[:-1]))]
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def get_valid_times(band, skypos, trange=None, detsize=1.1, verbose=0,
                    skyrange=None):
    """
    Given a sky position and (optional) extent, return all of the times
        periods containing spatially intersecting observations.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param trange: Minimum and maximum time (in GALEX time) to consider.

    :type trange: list

    :param detsize: The effective detector diameter, in degrees.

    :type detsize: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param skyrange: Values in degrees RA and Dec of a box around skypos that
        defines the extent of the region of interest.

    :type skyrange: list

    :returns: numpy.ndarray -- A sorted set of unique time stamps.
    """

    if not np.array(trange).tolist():
        trange = [1, 1000000000000]

    if len(np.shape(trange)) == 2:
        trange = trange[0]

    # Assembles sky positions on a grid within the targeted region.
    # [Future]: This is probably not an optimally efficient way to check an
    # entire region of sky for data, but it's not hugely dumb and does work...
    skypos_list = [skypos]
    if skyrange:
        # This massive construction with the hstack and separate calls to
        # linspace is to ensure that skypos (i.e. the target position) is
        # always uniquely searched.
        # In a perfect world, you would probably divice detsize by 2. The
        # detsize is divided by 3 to make sure to oversample the search just a
        # little bit.
        for r in np.unique(
                np.hstack([np.linspace(skypos[0],
                                       skypos[0]+skyrange[0]/2.,
                                       np.ceil(skyrange[0]/(detsize/3.)),
                                       endpoint=True),
                           np.linspace(skypos[0], skypos[0]-skyrange[0]/2.,
                                       np.ceil(skyrange[0]/(detsize/3.)),
                                       endpoint=True)])):
            for d in np.unique(
                    np.hstack([np.linspace(skypos[1],
                                           skypos[1]+skyrange[1]/2.,
                                           np.ceil(skyrange[1]/(detsize/3.)),
                                           endpoint=True),
                               np.linspace(skypos[1], skypos[1]-skyrange[1]/2.,
                                           np.ceil(skyrange[1]/(detsize/3.)),
                                           endpoint=True)])):
                skypos_list += [[r, d]]

    times = []
    for skypos in skypos_list:
        try:
            times = (list(times) +
                     list(np.array(gQuery.getArray(
                         gQuery.exposure_ranges(band, skypos[0], skypos[1],
                                                t0=trange[0], t1=trange[1],
                                                detsize=detsize),
                         verbose=verbose), dtype='float64')[:, 0]/tscale))
        except IndexError:
            if verbose:
                print("No exposure time available at {pos}".format(pos=skypos))
            return np.array([], dtype='float64')
        except TypeError:
            print("Is one of the inputs malformed?")
            raise
        except:
            raise

    # The following section should account for the fact that the use of band
    # in gQuery.exposure_ranges() doesn't seem to work like it should.
    # Specifically, it often returns time ranges as valid for FUV when there
    # are not actually any global events recorded by that detector.
    newtimes = []
    try:
        aspranges = distinct_tranges(np.sort(np.unique(times)))
    except IndexError:
        return np.array([], dtype='float64')
    except:
        raise

    for trange in aspranges:
        photontimes = (np.array(gQuery.getArray(
            gQuery.uniquetimes(band, trange[0], trange[1], flag=0)),
                                dtype='float64').flatten()/1000.)
        ix = np.where((photontimes.flatten() >= trange[0]) &
                      (photontimes.flatten() < trange[1]+1))
        if len(ix[0]):
            newtimes += np.arange(trange[0], trange[1]+1).tolist()

    return np.sort(np.unique(newtimes))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def fGetTimeRanges(band, skypos, trange=None, detsize=1.1, verbose=0,
                   maxgap=1., minexp=1., skyrange=None, maxgap_override=False):
    """
    Find the contiguous time ranges within a time range at a specific location.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param trange: Minimum and maximum time (in GALEX time) to consider.

    :type trange: list

    :param detsize: The effective detector diameter, in degrees.

    :type detsize: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param maxgap:  Maximum gap size, in seconds, for data to be considered
        contiguous.

    :type maxgap: float

    :param minexp: Minimum gap size, in seconds, for data to be considered
        contiguous.

    :type minexp: float

    :param skyrange: Values in degrees RA and Dec of a box around skypos that
        defines the extent of the region of interest.

    :type skyrange: list

    :param maxgap_override: Enables an experimental feature where maxgap
        can be less than one second.

    :type maxgap_override: bool

    :returns: numpy.ndarray -- A valid set of time ranges, accounting for
        minimum exposure lengths and maximum gaps.
	"""

    times = get_valid_times(band, skypos, trange=trange, detsize=detsize,
                            verbose=verbose, skyrange=skyrange)

    if not len(times):
        return np.array([[]], dtype='float64')

    if verbose:
        print_inline('Parsing ~'+str(len(times)-1)+' seconds of raw exposure.')

    # NOTE: The minimum meaningful maxgap is 1 second.
    if maxgap < 1 and not maxgap_override:
        raise ValueError('maxgap must be >=1 second')

    tranges = distinct_tranges(times, maxgap=maxgap)

    ix = np.where(np.array(tranges)[:, 1]-np.array(tranges)[:, 0] >= minexp)
    tranges = np.array(tranges)[ix].tolist()

    return np.array(tranges, dtype='float64')
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def stimcount_shuttered(band, trange, verbose=0, timestamplist=False):
    """
    Returns the stim count over a time range, excluding periods that the
        detector is considered shuttered (because of no non-NULL data).

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param trange: Minimum and maximum time (in GALEX time) to consider.

    :type trange: list

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param timestamplist: Global detector event timestamps.

    :type timestamplist: list

    :returns: int -- Total stim counts excluding shuttered time ranges.
    """

    try:
        t = (timestamplist if np.array(timestamplist).any() else
             np.array(gQuery.getArray(gQuery.uniquetimes(band, trange[0],
                                                         trange[1]),
                                      verbose=verbose),
                      dtype='float64')[:, 0]/gQuery.tscale)
    except IndexError: # Shutter this whole time range.
        if verbose:
            print('No data in {t0},{t1}'.format(t0=trange[0], t1=trange[1]))
        return 0

    times = np.sort(np.unique(np.append(t, trange)))
    tranges = distinct_tranges(times, maxgap=0.05)
    stimcount = 0

    for trange in tranges:
        stimcount += (gQuery.getValue(gQuery.stimcount(band, trange[0],
                                                       trange[1]),
                                      verbose=verbose) +
                      gQuery.getValue(gQuery.stimcount(band, trange[0],
                                                       trange[1], null=False),
                                      verbose=verbose))

    return stimcount
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def globalcount_shuttered(band, trange, verbose=0, timestamplist=False):
    """
    Global event counts over the time range, exluding shuttered periods (due to
        no non-NULL data).

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param trange: Minimum and maximum time (in GALEX time) to consider.

    :type trange: list

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param timestamplist: Global event time stamps.

    :type timestamplist: list

    :returns: int -- Total global counts excluding shuttered periods.
    """

    try:
        t = (timestamplist if np.array(timestamplist).any() else
             np.array(gQuery.getArray(gQuery.uniquetimes(band, trange[0],
                                                         trange[1], flag=True),
                                      verbose=verbose),
                      dtype='float64')[:, 0]/gQuery.tscale)
    except IndexError: # Shutter this whole time range.
        if verbose:
            print('No data in {t0},{t1}'.format(t0=trange[0], t1=trange[1]))
        return 0

    times = np.sort(np.unique(np.append(t, trange)))
    tranges = distinct_tranges(times, maxgap=0.05)
    nonnullevents, nullevents = 0, 0

    for trange in tranges:
        nullevents += gQuery.getValue(
            gQuery.deadtime2(band, trange[0], trange[1]), verbose=verbose)

        nonnullevents += gQuery.getValue(gQuery.deadtime1(band, trange[0],
                                                          trange[1]),
                                         verbose=verbose)

    return nullevents+nonnullevents
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def compute_shutter(band, trange, verbose=0, shutgap=0.05,
                    timestamplist=False):
    """
    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param trange: Minimum and maximum time (in GALEX time) to consider.

    :type trange: list

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param shutgap: Amount of time, in seconds, that defines the minimum gap in
        observation that corresponds to a 'shutter' (not a true exposure time).

    :type shutgap: float

    :param timestamplist: Global event time stamps.

    :type timestamplist: list

    :returns: numpy.ndarray -- The total shutter time, in seconds, during the
        specified time range.
    """

    try:
        t = (timestamplist if np.array(timestamplist).any() else
             np.array(gQuery.getArray(gQuery.uniquetimes(band, trange[0],
                                                         trange[1], flag=True),
                                      verbose=verbose),
                      dtype='float64')[:, 0]/gQuery.tscale)
    except IndexError: # Shutter this whole time range if there's no data
        return trange[1]-trange[0]

    t = np.sort(np.unique(np.append(t, trange)))
    ix = np.where(t[1:]-t[:-1] >= shutgap)

    return np.array(t[1:]-t[:-1])[ix].sum()
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def empirical_deadtime(band, trange, verbose=0, feeclkratio=0.966,
                       timestamplist=False):
    """
    Calculate empirical deadtime (per global count rate) using revised
        formulas. Restricts integration of global counts to non-shuttered time
        periods.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param trange: Minimum and maximum time (in GALEX time) to consider.

    :type trange: list

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param feeclkratio: A scaling parameter for the Front End Electronics clock.

    :type feeclkratio: float

    :param timestamplist: Global even time stamps.

    :type timestamplist: list

    :returns: float -- The empirical deadtime ratio.
    """

    model = {'NUV':[-0.000434730599193, 77.217817988],
             'FUV':[-0.000408075976406, 76.3000943221]}

    rawexpt = trange[1]-trange[0]-compute_shutter(band, trange,
                                                  timestamplist=timestamplist)
    gcr = globalcount_shuttered(band, trange,
                                timestamplist=timestamplist)/rawexpt

    refrate = model[band][1]/feeclkratio

    scr = model[band][0]*gcr+model[band][1]

    return 1 - scr/feeclkratio/refrate
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def exposure(band, trange, verbose=0):
    """
    Calculate the effective exposure time in a period, in seconds, accounting
        for shutter and deadtime. Does not account for actual sky coverage of
        the telescope during the time period queried (see: compute_exptime()
        below).

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param trange: Minimum and maximum time (in GALEX time) to consider.

    :type trange: list

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :returns: float -- The effective exposure time, in seconds.
    """

    rawexpt = trange[1]-trange[0]

    if rawexpt == 0.:
        return 0.

    try:
        t = (np.array(gQuery.getArray(
            gQuery.uniquetimes(band, trange[0], trange[1], flag=True),
            verbose=verbose), dtype='float64')[:, 0]/gQuery.tscale)
    except IndexError: # Shutter this whole time range.
        if verbose:
            print('No data in {t0},{t1}'.format(t0=trange[0], t1=trange[1]))
        return 0.

    shutter = compute_shutter(band, trange, verbose=verbose, timestamplist=t)

    deadtime = empirical_deadtime(band, trange, verbose=verbose, timestamplist=t)

    return (rawexpt-shutter)*(1.-deadtime)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def compute_exptime(band, tr, verbose=0, skypos=None, detsize=1.25,
                    coadd=False):
    """
    Compute the total effective exposure time, in seconds, accounting for
        shutter and deadtime _and_ detector size (i.e. effective coverage).

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param tr: Pairs of minimum and maximum times (in GALEX time) to consider.

    :type tr: list

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param skypos: The right ascension and declination of interest, in degrees.

    :type skypos: list

    :param detsize: The effective detector diameter, in degrees.

    :type detsize: float

    :param coadd: Should the effective exposure time be calculated across all
        time ranges, e.g., a coadded effective exposure time.

    :type coadd: bool

    :returns: list
    """

    if len(np.shape(tr)) == 1:
        tr = [tr]

    if not len(np.shape(tr)) == 2:
        raise ValueError('tr array has {d} dimensions,'
                         ' must be <=2'.format(d=len(np.shape(tr))))

    if skypos:
        exptime = []
        for trange in tr:
            tranges = fGetTimeRanges(band, skypos, verbose=verbose,
                                     trange=trange, detsize=detsize).tolist()
            if np.array(tranges).any():
                exptime += [sum(
                    exposure(band, trange, verbose=verbose)
                    for trange in tranges)]
            else:
                exptime += [0.]
    else:
        exptime = [exposure(band, trange, verbose=verbose)
                   for trange in tr]

    return (([sum(exptime)] if coadd else exptime)
            if np.array(tr).any() else [0.])
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def get_mcat_data(skypos, rad):
    """
    Return visit-level MCAT sources and their catalog values within a give
        radius of the specified sky position.

    :param skypos: The right ascension and declination, in degrees, around
        which to search for MCAT sources.

    :type skypos: list

    :param rad: The radius within which to search for MCAT sources,
        in degrees.

    :type rad: float

    :returns: dict -- The MCAT sources within the specified radius.
    """

    # Try once with the default radius.
    out = np.array(gQuery.getArray(gQuery.mcat_visit_sources(skypos[0],
                                                             skypos[1], rad)))
    # If no MCAT sources found, try again with a radius 5 times bigger.
    if len(out) == 0:
        out = np.array(gQuery.getArray(gQuery.mcat_visit_sources(skypos[0],
                                                                 skypos[1],
                                                                 rad*5.)))
    # [Future]: The APER entries should really be generated.
    try:
        return {'objid':np.array(out[:, 0], dtype='int64'),
                'ra':np.array(out[:, 1], dtype='float32'),
                'dec':np.array(out[:, 2], dtype='float32'),
                'NUV':{'mag':np.array(out[:, 3], dtype='float32'),
                       't0':np.array(out[:, 46], dtype='float64')-GPSSECS,
                       't1':np.array(out[:, 47], dtype='float64')-GPSSECS,
                       'skybg':np.array(out[:, 6], dtype='float32'),
                       'expt':np.array(out[:, 11], dtype='float32'),
                       'fwhm':np.array(out[:, 8], dtype='float32'),
                       'obssecs':np.array(out[:, 40], dtype='float64'),
                       'artifact':np.array(out[:, 42], dtype='int32'),
                       'ra':np.array(out[:, 50], dtype='float32'),
                       'dec':np.array(out[:, 51], dtype='float32'),
                       1:{'mag':np.array(out[:, 19],
                                         dtype='float32')+zpmag('NUV'),
                          'err':np.array(out[:, 33], dtype='float32')},
                       2:{'mag':np.array(out[:, 20],
                                         dtype='float32')+zpmag('NUV'),
                          'err':np.array(out[:, 34], dtype='float32')},
                       3:{'mag':np.array(out[:, 21],
                                         dtype='float32')+zpmag('NUV'),
                          'err':np.array(out[:, 35], dtype='float32')},
                       4:{'mag':np.array(out[:, 22],
                                         dtype='float32')+zpmag('NUV'),
                          'err':np.array(out[:, 36], dtype='float32')},
                       5:{'mag':np.array(out[:, 23],
                                         dtype='float32')+zpmag('NUV'),
                          'err':np.array(out[:, 37], dtype='float32')},
                       6:{'mag':np.array(out[:, 24],
                                         dtype='float32')+zpmag('NUV'),
                          'err':np.array(out[:, 38], dtype='float32')},
                       7:{'mag':np.array(out[:, 25],
                                         dtype='float32')+zpmag('NUV'),
                          'err':np.array(out[:, 39], dtype='float32')}},

                'FUV':{'mag':np.array(out[:, 4], dtype='float32'),
                       't0':np.array(out[:, 44], dtype='float64')-GPSSECS,
                       't1':np.array(out[:, 45], dtype='float64')-GPSSECS,
                       'skybg':np.array(out[:, 7], dtype='float32'),
                       'expt':np.array(out[:, 10], dtype='float32'),
                       'fwhm':np.array(out[:, 9], dtype='float32'),
                       'obssecs':np.array(out[:, 41], dtype='float64'),
                       'artifact':np.array(out[:, 43], dtype='int32'),
                       'ra':np.array(out[:, 48], dtype='float32'),
                       'dec':np.array(out[:, 49], dtype='float32'),
                       1:{'mag':np.array(out[:, 12],
                                         dtype='float32')+zpmag('FUV'),
                          'err':np.array(out[:, 26], dtype='float32')},
                       2:{'mag':np.array(out[:, 13],
                                         dtype='float32')+zpmag('FUV'),
                          'err':np.array(out[:, 27], dtype='float32')},
                       3:{'mag':np.array(out[:, 14],
                                         dtype='float32')+zpmag('FUV'),
                          'err':np.array(out[:, 28], dtype='float32')},
                       4:{'mag':np.array(out[:, 15],
                                         dtype='float32')+zpmag('FUV'),
                          'err':np.array(out[:, 29], dtype='float32')},
                       5:{'mag':np.array(out[:, 16],
                                         dtype='float32')+zpmag('FUV'),
                          'err':np.array(out[:, 30], dtype='float32')},
                       6:{'mag':np.array(out[:, 17],
                                         dtype='float32')+zpmag('FUV'),
                          'err':np.array(out[:, 31], dtype='float32')},
                       7:{'mag':np.array(out[:, 18],
                                         dtype='float32')+zpmag('FUV'),
                          'err':np.array(out[:, 32], dtype='float32')}}}
    except IndexError:
        # If there are STILL no detections, then pass a dict with empty values.
        # A default set of values will then be used.
        return {'objid':None, 'ra':None, 'dec':None, 'NUV':None, 'FUV':None}
    except:
        raise
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def exp_from_objid(objid):
    """
    Return the effective exposure time, start time, and end time from the MCAT
        for a given GALEX object ID.

    :param objid: GALEX object ID.

    :type objid: int

    :returns: dict -- The FUV and NUV effective exposure time and start/stop
        times.
    """

    out = np.array(gQuery.getArray(gQuery.mcat_objid_search(objid)))

    return {'NUV':
                {'expt':np.array(out[:, 7], dtype='float')[0],
                 't0':np.array(out[:, 9], dtype='float64')[0]-GPSSECS,
                 't1':np.array(out[:, 10], dtype='float64')[0]-GPSSECS},
            'FUV':{'expt':np.array(out[:, 8], dtype='float')[0],
                   't0':np.array(out[:, 11], dtype='float64')[0]-GPSSECS,
                   't1':np.array(out[:, 12], dtype='float64')[0]-GPSSECS}}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def obstype(t, obsdata=None):
    """
    Determines the type of dither pattern.
    """
    if not obsdata:
        obsdata = gQuery.getArray(gQuery.obstype_from_t(t))
    try:
        return str(obsdata[0][0])
    except IndexError:
        return "Unknown"
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def legnum(t, obsdata=None):
    """
    Returns the lef number.
    """
    if not obsdata:
        obsdata = gQuery.getArray(gQuery.obstype_from_t(t))
    try:
        return obsdata[0][5]
    except IndexError:
        return "Unknown"
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def obstype_from_objid(objid):
    """
    Return the number of legs and petal value for a given GALEX object ID.

    :param objid: GALEX object ID.

    :type objid: int

    :returns: tuple -- A two-element tuple containing the number of legs and
        the petal value, which can be used to infer the observation
        type/strategy.
    """

    out = gQuery.getArray(gQuery.obstype(objid))

    (survey, tilename, photoextractid, petal, nlegs, leg, eclipse, img,
     subvis) = np.array(out)[0]

    return int(nlegs), int(petal)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def mcat_skybg(band, skypos, radius, verbose=0, trange=None, mcat=None,
               searchradius=0.1):
    """
    Estimate the sky background using the MCAT 'skybg' for nearby sources.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param radius: The radius in which to search for MCAT sources in degrees.

    :type radius: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param trange: Minimum and maximum time (in GALEX time) to consider.

    :type trange: list

    :returns: float -- The estimated sky background in the photometric
        aperture, in counts per second.
    """
    # Search the visit-level MCAT for nearby detections.
    # Unless the MCAT data has already been handed off for detection purposes.
    if not mcat:
        mcat = get_mcat_data(skypos, searchradius)
    try:
        # Find the distance to each source.
        dist = np.array([angularSeparation(skypos[0], skypos[1], a[0], a[1])
                         for a in zip(mcat[band]['ra'], mcat[band]['dec'])])
    except TypeError:
        print_inline(
            'No {b} MCAT sources within {r} degrees of {p}'.format(
                b=band, r=searchradius, p=skypos))
        return np.nan

    # Find visits that overlap in time.
    if not trange:
        tix = (np.array(list(range(len(mcat[band]['mag']))), dtype='int32'),)
    else:
        tix = np.where(
            ((trange[0] >= mcat[band]['t0']) & (trange[0] <= mcat[band]['t1'])) |
            ((trange[1] >= mcat[band]['t0']) & (trange[1] <= mcat[band]['t1'])) |
            ((trange[0] <= mcat[band]['t0']) & (trange[1] >= mcat[band]['t1'])))

    if not len(tix[0]):
        print_inline('No concurrent {b} MCAT source nearby.'.format(b=band))
        return np.nan # Might not be the preferred behavior here.

    ix = np.where(dist[tix] == min(dist[tix]))

    skybg = mcat[band]['skybg'][tix][ix]

    # This should rarely happen, but sometimes there's a duplicate entry in
    # the visit-level MCAT.
    if len(skybg) > 1:
        # If the skybg array is all the same value, it's a duplicate.
        if np.all(skybg == skybg[0]):
            skybg = np.asarray([skybg[0]])
        else:
            skybg = np.asarray([np.median(skybg)])

    return skybg[0]*area(radius*60.*60.)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def get_mags(band, ra0, dec0, radius, maglimit, mode='coadd',
             zpmag={'NUV':20.08, 'FUV':18.82}, verbose=0):
    """
    Given RA, Dec and search radius, searches the coadd MCAT for sources.
        Returns a dict() which contains magnitudes for all of the APER settings.
        Note: Visit mode returns a lot more sources, more slowly than coadd mode
        given the same search parameters. You should probably use smaller search
        radii in visit mode. If you're just trying to find unique sources in a
        large region, use coadd mode and then pass the result through the
        parse_unique_sources() function contained in this module.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: The right ascension, in degrees, around which to search.

    :type ra0: float

    :param dec0: The declination, in degrees, around which to search.

    :type dec0: float

    :param radius: The size of the search radius for MCAT sources in degrees.

    :type radius: float

    :param maglimit: The NUV faint limit to return MCAT sources for.

    :type maglimit: float

    :param mode: Specify whether to return MCAT sources from the 'visit' or
        'coadd' catalog.

    :type mode: str

    :param zpmag: The zero-point magnitude in FUV and NUV.

    :type zpmag: dict

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :returns: dict -- The set of magnitudes from different apertures for sources
        in the MCAT, centered around the specified coordinate.
    """

    zpf, zpn = zpmag['FUV'], zpmag['NUV']

    if mode == 'coadd':
        out = np.array(gQuery.getArray(gQuery.mcat_sources(band, ra0, dec0,
                                                           radius,
                                                           maglimit=maglimit),
                                       verbose=verbose))
        if not len(out):
            print("Warning: No sources found!")
            return None
        return {'ra':out[:, 0], 'dec':out[:, 1],
                'FUV':{'mag':out[:, 3], 1:out[:, 9]+zpf, 2:out[:, 10]+zpf,
                       3:out[:, 11]+zpf, 4:out[:, 12]+zpf, 5:out[:, 13]+zpf,
                       6:out[:, 14], 7:out[:, 15]+zpf},
                'NUV':{'mag':out[:, 2], 1:out[:, 16]+zpn, 2:out[:, 17]+zpn,
                       3:out[:, 18]+zpn, 4:out[:, 19]+zpn, 5:out[:, 20]+zpn,
                       6:out[:, 21]+zpn, 7:out[:, 22]+zpn}}
    elif mode == 'visit':
        return get_mcat_data([ra0, dec0], radius)
    else:
        print("mode must be in [coadd,visit]")
        return None
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def find_nearest_mcat(band, skypos, radius, maglimit=30.):
    """
    Given a sky position and a search radius, find the nearest MCAT source
        and return its position and magnitude in specified band.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: Two element array of RA and Dec in decimal degrees.

    :type skypos: array

    :param radius: Search radius in decimal degrees.

    :type radius: float

    :param maglimit: The NUV faint limit to return MCAT sources for.

    :type maglimit: float
    """

    data = get_mags(band, skypos[0], skypos[1], radius, 30)
    if not data:
        return None

    separation = [angularSeparation(skypos[0], skypos[1], a[0], a[1])
                  for a in zip(data['ra'], data['dec'])]
    minsep = np.where(separation == min(separation))

    return {'mag':data[band]['mag'][minsep][0],
            'skypos':np.array(list(zip(data['ra'],
                                       data['dec'])))[minsep][0].tolist(),
            'distance':min(separation)}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def parse_unique_sources(ras, decs, margin=0.001):
    """
    Iteratively returns unique sources based upon a margin within
        which two sources should be considered the same sources. Is a little
        bit sensitive to the first entry and could probably be written to be
        more robust, but works well enough.

    :param ras: Set of right ascensions, in degrees.

    :type ras: numpy.ndarray

    :param decs: Set of declinations, in degrees.

    :type decs: numpy.ndarray

    :param margin: Buffer size when determining matches, in degrees.

    :type margin: float

    :returns: list -- The indices corresponding to unique sources.
    """

    skypos = list(zip(ras, decs))

    for i, pos in enumerate(skypos):
        ix = np.where(angularSeparation(pos[0], pos[1], ras, decs) <= margin)
        skypos[i] = [ras[ix].mean(), decs[ix].mean()]
        a = skypos

    b = []

    for i in a:
        if not i in b:
            b += [i]

    return b
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def find_unique_sources(band, ra0, dec0, searchradius, maglimit=20.0,
                        margin=0.001, verbose=0):
    """
    Locates nominally unique (via crossmatch) GALEX sources in the MCAT
        near a sky position of interest.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: The right ascension, in degrees, around which to search.

    :type ra0: float

    :param dec0: The declination, in degrees, around which to search.

    :type dec0: float

    :param searchradius: The size of the radius to search for unique sources,
        in degrees.

    :type searchradius: float

    :param maglimit: The NUV faint limit to return unique sources for.

    :type maglimit: float

    :param margin: Buffer size when determining matches, in degrees.

    :type margin: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :returns: numpy.ndarray -- The set of unique sources.
    """

    coadds = get_mags(band, ra0, dec0, searchradius, maglimit, mode='coadd',
                      verbose=verbose)

    if not coadds:
        return None
    else:
        return np.array(parse_unique_sources(coadds['ra'], coadds['dec'],
                                             margin=margin))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def avg_sources(band, skypos, radius=0.001, maglimit=20.0, verbose=0,
                catalog='MCAT'):
    """
    Return the mean position of sources within the search radius.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param radius: The radius within which to search for MCAT sources,
        in degrees?

    :type radius: float

    :param maglimit: The NUV faint limit to return unique sources for.

    :type maglimit: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param catalog: The name of the catalog to query.

    :type catalog: str

    :returns: tuple -- A three-element tuple containing the mean RA, mean DEC,
        and mean FWHM of sources within the search radius.
    """

    out = np.array(gQuery.getArray(gQuery.mcat_sources(band, skypos[0],
                                                       skypos[1], radius,
                                                       maglimit=maglimit),
                                   verbose=verbose))

    ix = np.where(out[:, -2] > 0) if band == 'NUV' else np.where(out[:, -1] > 0)

    fwhm = out[ix, -2].mean() if band == 'NUV' else out[ix, -1].mean()

    return out[ix, 0].mean(), out[ix, 1].mean(), round(fwhm, 4)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def nearest_source(band, skypos, radius=0.01, maglimit=20.0, verbose=0,
                   catalog='MCAT'):
    """
    Return targeting parameters for the nearest MCAT source to a position.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param radius: The radius within which to search for the nearest MCAT
        source, in degrees.

    :type radius: float

    :param maglimit: The NUV faint limit to return unique sources for.

    :type maglimit: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param catalog: The name of the catalog to query.

    :type catalog: str

    :returns: tuple -- A three-element tuple containing the mean RA, mean DEC,
        and mean FWHM of the nearest sources within the search radius.
    """

    out = np.array(gQuery.getArray(gQuery.mcat_sources(band, skypos[0],
                                                       skypos[1], radius,
                                                       maglimit=maglimit),
                                   verbose=verbose))

    if not len(out) and band == 'FUV':
        if verbose:
            print("No nearby MCAT source found in FUV. Trying NUV...")
        band = 'NUV'
        out = np.array(gQuery.getArray(gQuery.mcat_sources(band, skypos[0],
                                                           skypos[1], radius,
                                                           maglimit=maglimit),
                                       verbose=verbose))

    if not len(out) and band == 'NUV':
        if verbose:
            print("No nearby MCAT source found. Using input sky position.")
        return skypos[0], skypos[1], 0.01

    dist = angularSeparation(out[:, 0], out[:, 1], skypos[0], skypos[1])

    if verbose > 1:
        print("Finding nearest among "+str(len(dist))+" nearby sources.")

    # Note that this doesn't cope with multiple entries for the same source.
    s = out[np.where(dist == dist.min())][0]

    return avg_sources(band, [s[0], s[1]], verbose=verbose)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def nearest_distinct_source(band, skypos, radius=0.1, maglimit=20.0, verbose=0,
                            catalog='MCAT'):
    """
    Return parameters for the nearest non-targeted source.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param radius: The radius within which to search for the nearest MCAT
        source, in degrees.

    :type radius: float

    :param maglimit: The NUV faint limit to return unique sources for.

    :type maglimit: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param catalog: The name of the catalog to query.

    :type catalog: str

    :returns: numpy.ndarray -- Catalog values for the nearest non-targeted
        source.
    """

    out = np.array(gQuery.getArray(gQuery.mcat_sources(band, skypos[0],
                                                       skypos[1], radius,
                                                       maglimit=maglimit),
                                   verbose=verbose))

    dist = angularSeparation(out[:, 0], out[:, 1], skypos[0], skypos[1])

    ix = np.where(dist > 0.005)

    return np.array(out)[ix][np.where(dist[ix] == dist[ix].min())][0]
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def suggest_bg_radius(band, skypos, radius=0.1, maglimit=20.0, verbose=0,
                      catalog='MCAT'):
    """
    Returns a recommended background radius based upon the positions and FWHM of
        nearby sources in the MCAT.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param radius: The radius within which to search for MCAT  sources,
        in degrees?

    :type radius: float

    :param maglimit: The NUV faint limit to return unique sources for.

    :type maglimit: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param catalog: The name of the catalog to query.

    :type catalog: str

    :returns: float -- The suggested radius to define the background with.
    """

    nearest = nearest_distinct_source(band, skypos, verbose=verbose)

    dist = angularSeparation(nearest[0], nearest[1], skypos[0], skypos[1])

    return round(dist-3*nearest[-2 if band == 'NUV' else -1], 4)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def optimize_annulus(optrad, outann, verbose=0):
    """
    Suggest optiumum annulus dimensions.

    :param optrad: The optimal photometric aperture radius, in degrees.

    :type optrad: float

    :param outann: The outer annulus to test, in degrees.

    :type outann: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :returns: float - The suggested outer annulus for background calculation.
    """

    if outann <= round(2*optrad, 4):
        print("Warning: There are known sources within the background annulus.")
        print("Use --hrbg to mask these out. (Will increase run times.)")

    return round(1.2*optrad, 4), round(2*optrad, 4)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def suggest_parameters(band, skypos, verbose=0):
    """
    Provide suggested coordinates and photometric apertures for a source
        given the location of known MCAT sources nearby.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :returns: tuple -- A five-element tuple containing the suggested right
        ascension, declination, photometric aperture, inner annulus, and outer
        annulus, all in degrees.
    """

    mcat = get_mcat_data(skypos, 0.01)
    ix = np.where((mcat[band]['mag'] > 0) & (mcat[band]['fwhm'] > 0))
    pos, fwhm = None, None

    if mcat['objid'].any(): # There is a known star at the target position!
        pos = [mcat['ra'][ix].mean(), mcat['dec'][ix].mean()]
        radius = 2*mcat[band]['fwhm'][ix].mean()
        if verbose:
            print('Recentering on {pos}.'.format(pos=pos))
            print('Using aperture radius of {rad} degrees.'.format(rad=fwhm))
    else: # There is no known star at the target position...
        pos = skypos
        radius = aper2deg(4)
    annulus = [3*radius, 5*radius]

    return pos[0], pos[1], radius, annulus[0], annulus[1]
# ------------------------------------------------------------------------------
