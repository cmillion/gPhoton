"""
.. module:: curvetools

   :synopsis: Methods for creating light curves.
   @CHASE - elaborate on this please.@

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

import numpy as np
import pandas as pd
import os
# gPhoton specific
import gQuery
from gQuery import tscale
import MCUtils as mc
import dbasetools as dbt
import galextools as gxt
# @CHASE - Looks like there's a bunch of methods in cal/__init.py__, shouldn't
# these be put into a module?@
import cal

# ------------------------------------------------------------------------------
def gphot_params(band, skypos, radius, annulus=None, verbose=0, detsize=1.25,
                 stepsz=None, trange=None):
    """
    Populate a dict() with parameters that are constant over all bins.

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: A two-element list containing the RA and DEC.

    :type skypos: list

    :param radius: The photometric aperture, in degrees. @CHASE - confirm units@

    :type radius: float

    :param annulus: A two-element list containing the inner and outer radius
    to use for background subtraction during aperture photometry, in degrees.
    @CHASE - confirm units.@

    :type annulus: list

    :param verbose: Level of verbosity, 0 = minimum verbosity.

    :type verbose: int @CHASE - The default was float, I made this an int.@

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :param stepsz: Size of time bin to use, in seconds.

    :type stepsz: float

    :param trange: The start and end timestamps to consider. @CHASE - JD, GALEX
    time?@

    :type trange: list

    :returns: dict -- The set of parameters that are constant across all bins.
    """

    return {'band':band, 'ra0':skypos[0], 'dec0':skypos[1], 'skypos':skypos,
            'trange':trange, 'radius':radius, 'annulus':annulus,
            'stepsz':stepsz, 'verbose':verbose, 'detsize':detsize,
            'apcorrect1':gxt.apcorrect1(radius, band),
            'apcorrect2':gxt.apcorrect2(radius, band),
            'detbg':gxt.detbg(mc.area(radius), band)}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def xieta2colrow(xi, eta, band, detsize=1.25):
    """
    Convert detector xi and eta into col and row.

    :param xi: @CHASE - What is 'xi'?@

    :type xi: @CHASE - float/list/numpy.ndarray?@

    :param eta: @CHASE - What is 'eta'?@

    :type eta: @CHASE - float/list/numpy.ndarray?@

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

    # @CHASE - Should this note be moved out of the source code for clarity?@
    # You could theoretically drop a cut on detector position / detsize here...
    # Also, is this cut absolutely necessary? I think it's already been taken
    # care of by the flag==0 assertion in the SQL query.
    # cut = ((col > 0.) & (col < flat.shape[0]-1) &
    #       (row > 0.) & (row < flat.shape[1]-1))
    # cut = np.where(ix == True)
    # ix = np.where((1.25/800.)*mc.distance(col,row,400,400)=detsize)
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

    :param ra0: Right ascension, in degrees, of the center of the field-of-view.
    @CHASE - confirm please.@

    :type ra0: float

    :param dec0: Declination, in degrees, of the center of the field-of-view.
    @CHASE - confirm please.@

    :type dec0: float

    :param tranges: Set of time ranges to retrieve the photon events.
    @CHASE - in GALEX time?@

    :type tranges: list

    :param radius: The photometric aperture, in degrees. @CHASE - Confirm
    this is in degrees. Also, the retrieval only gets those that are within
    the radius and not the outer annulus?@

    :type radius: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param colnames: Labels of the columns found in the photon event file.

    :type colnames: list

    :returns: dict -- The set of photon events and their properties.
    """

    # [Future]: Consider moving this method to 'dbasetools'.

    if verbose:
        print 'Reading photon list file: {f}'.format(f=photonfile)

    data = pd.io.parsers.read_csv(photonfile, names=colnames)
    ra, dec = np.array(data['ra']), np.array(data['dec'])
    angsep = mc.angularSeparation(ra0, dec0, ra, dec)

    ix = np.array([])
    for trange in tranges:
        if verbose:
            print trange
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
def query_photons(band, ra0, dec0, tranges, radius, verbose=0, flag=0):
    """
    Retrieve photons within an aperture from the database.

    :param band: Name of the band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: Right ascension, in degrees, of the center of the field-of-view.
    @CHASE - confirm please.@

    :type ra0: float

    :param dec0: Declination, in degrees, of the center of the field-of-view.
    @CHASE - confirm please.@

    :type dec0: float

    :param tranges: Set of time ranges to retrieve the photon events.
    @CHASE - in GALEX time?@

    :type tranges: list

    :param radius: The photometric aperture, in degrees. @CHASE - Confirm
    this is in degrees. Also, the retrieval only gets those that are within
    the radius and not the outer annulus?@

    :type radius: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param flag: @CHASE - What is 'flag'?@

    :type flag: int

    :returns: dict -- The set of photon events with their properties.
    """

    # [Future]: This should be moved to 'dbasetools'.
    stream = []
    if verbose:
        print "Retrieving photons within {rad} degrees of [{r}, {d}]".format(
            rad=radius, r=ra0, d=dec0)
    for trange in tranges:
        if verbose:
            mc.print_inline(" and between {t0} and {t1}.".format(t0=trange[0],
                                                                 t1=trange[1]))
        thisstream = gQuery.getArray(
            gQuery.allphotons(band, ra0, dec0, trange[0], trange[1], radius,
                              flag=flag), verbose=verbose, retries=100)
        stream.extend(thisstream)

    stream = np.array(stream, 'f8').T
    colnames = ['t', 'ra', 'dec', 'xi', 'eta', 'x', 'y']
    dtypes = ['f8', 'f8', 'f8', 'f4', 'f4', 'f4', 'f4']
    cols = map(np.asarray, stream, dtypes)

    events = dict(zip(colnames, cols))

    # Adjust the timestamp by tscale.
    events['t'] /= tscale

    return events
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def pullphotons(band, ra0, dec0, tranges, radius, events={}, verbose=0,
                photonfile=None, flag=0):
    """
    @CHASE - Please provide descriptor.@

    :param band: Name of the band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: Right ascension, in degrees, of the center of the field-of-view.
    @CHASE - confirm please.@

    :type ra0: float

    :param dec0: Declination, in degrees, of the center of the field-of-view.
    @CHASE - confirm please.@

    :type dec0: float

    :param tranges: Set of time ranges to retrieve the photon events.
    @CHASE - In GALEX time?@

    :type tranges: list

    :param radius: The photometric aperture, in degrees. @CHASE - Confirm
    this is in degrees. Also, the retrieval only gets those that are within
    the radius and not the outer annulus?@

    :type radius: float

    :param events: Set of photon events. @CHASE - Why is this passed as an
    argument? It's defined within the method and returned.@

    :type events: dict

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param photonfile: Name of photon event file to use (if reading from disk).

    :type photonfile: str

    :param flag: @CHASE - What is 'flag'?@

    :type flag: int

    :returns: dict -- The set of photon events and their properties.
    """

    if photonfile:
        events = read_photons(photonfile, ra0, dec0, tranges, radius,
                              verbose=verbose)
    else:
        events = query_photons(band, ra0, dec0, tranges, radius,
                               verbose=verbose, flag=flag)

    events = hashresponse(band, events, verbose=verbose)

    return events
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def aperture_error(counts, expt, bgcounts=0):
    """
    @CHASE - Please provide descriptor.@

    :param counts: Total counts within the aperture. @CHASE - confirm please.@

    :type counts: int @CHASE - Assume not a float?@

    :param expt: The exposure time in seconds.

    :type extp: float @CHASE - Is this a float or int?@

    :param bgcounts: The total background counts within the aperture.

    :type bgcounts: int @CHASE - Assume not a float?@

    :returns: float -- The error in the counts within the aperture.
    """

    return np.sqrt(counts/expt**2+bgcounts/expt**2)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def bg_sources(band, ra0, dec0, radius, margin=0.001):
    """
    @CHASE - Please provide descriptor.@

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: Right ascension, in degrees, of the center of the field-of-view.
    @CHASE - confirm please.@

    :type ra0: float

    :param dec0: Declination, in degrees, of the center of the field-of-view.
    @CHASE - confirm please.@

    :type dec0: float

    :param radius: The search radius to find MCAT sources within, in degrees.
    @CHASE - please confirm description and units.@

    :type radius: float

    :param margin: Extra search margin when looking for MCAT sources, in
    degrees. @CHASE - confirm description and units.@

    :type margin: float

    :returns: dict -- The RA, DEC, Full-Width-Half-Maximum, Mag. Limit, and
    search radius. @CHASE - please confirm what 'radius' is in return dict.@
    """

    sources = gQuery.getArray(gQuery.mcat_sources(band, ra0, dec0,
                                                  radius+margin,
                                                  maglimit=maskdepth))

    try:
        return {'ra':np.float32(np.array(sources)[:, 0]),
                'dec':np.float32(np.array(sources)[:, 1]),
                'fwhm':np.float32(np.array(sources)[:, 7:9]),
                'radius':radius}
    except IndexError:
        return {'ra':np.array([]), 'dec':np.array([]), 'fwhm':np.array([]),
                'maglimit':maskdepth, 'radius':radius}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def bg_mask_annulus(band, ra0, dec0, annulus, ras, decs, responses):
    """
    @CHASE - Please provide descriptor.@

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: Right ascension, in degrees, of the center of the aperture.
    @CHASE - confirm please.@

    :type ra0: float

    :param dec0: Declination, in degrees, of the center of the aperture.
    @CHASE - confirm please.@

    :type dec0: float

    :param annulus: The inner and outer radius of the annulus used to define
    the background, in degrees.

    :type annulus: list

    :param ras: Right ascension of sources to check whether they lie in the
    annulus, in degrees.

    :type ras: numpy.ndarray @CHASE - confirm data type is not list.@

    :param decs: Declination of sources to check whether they lie in the
    annulus, in degrees.

    :type decs: numpy.ndarray @CHASE - confirm data type is not list.@

    :param responses: Set of response values for the sources that will be
    checked whether they lie in the annulus.

    :type responses: numpy.ndarray @CHASE - confirm data type is not list.@

    :returns: tuple -- A three-element tuple containing the RAs, DECs, and
    responses for sources that lie within the specified annulus.
    """

    # @CHASE - What happens if np.where() returns no matches? Should be caught?@

    ix = np.where((mc.angularSeparation(ra0, dec0, ras, decs) >= annulus[0]) &
                  (mc.angularSeparation(ra0, dec0, ras, decs) <= annulus[1]))

    return ras[ix], decs[ix], responses[ix]
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def bg_mask_sources(band, ra0, dec0, ras, decs, responses, sources,
                    maskradius=1.5):
    """
    @CHASE - Please provide descriptor.@

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: Right ascension, in degrees, of the center of the aperture.
    @CHASE - confirm please.@

    :type ra0: float

    :param dec0: Declination, in degrees, of the center of the aperture.
    @CHASE - confirm please.@

    :type dec0: float

    :param ras: Right ascension of sources to check whether they lie in the
    annulus, in degrees.

    :type ras: numpy.ndarray @CHASE - confirm data type is not list.@

    :param decs: Declination of sources to check whether they lie in the
    annulus, in degrees.

    :type decs: numpy.ndarray @CHASE - confirm data type is not list.@

    :param responses: Set of response values for the sources that will be
    checked whether they lie in the annulus.

    :type responses: numpy.ndarray @CHASE - confirm data type is not list.@

    :param sources: The set of sources to check. @CHASE - please update.@

    :type sources: dict @CHASE - please confirm data type.@

    :param maskradius: The value to use for the mask radius.

    :type maskradius: float

    :returns: tuple -- A three-element tuple containing the RAs, DECs, and
    responses for sources beyond the mask radius. @CHASE - Please update
    this description.@
    """

    # At present, masks to 1.5 sigma where FWHM = 2.3548*sigma.
    for i in range(len(sources['ra'])):
        ix = np.where(
            mc.angularSeparation(
                sources['ra'][i], sources['dec'][i], ras, decs) >=
            (maskradius/2.3548)*np.median(sources['fwhm'][i, :]))
        ras, decs, responses = ras[ix], decs[ix], responses[ix]

    return ras, decs, responses
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def bg_mask(band, ra0, dec0, annulus, ras, decs, responses, sources):
    """
    @CHASE - Please provide descriptor.@

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: Right ascension, in degrees, of the center of the aperture.
    @CHASE - confirm please.@

    :type ra0: float

    :param dec0: Declination, in degrees, of the center of the aperture.
    @CHASE - confirm please.@

    :type dec0: float

    :param annulus: Size of the inner and outer annuli, in degrees. @CHASE -
    verify unit please.@

    :type annulus: float

    :param ras: Right ascension of sources to check whether they lie in the
    annulus, in degrees.

    :type ras: numpy.ndarray @CHASE - confirm data type is not list.@

    :param decs: Declination of sources to check whether they lie in the
    annulus, in degrees.

    :type decs: numpy.ndarray @CHASE - confirm data type is not list.@

    :param responses: Set of response values for the sources that will be
    checked whether they lie in the annulus.

    :type responses: numpy.ndarray @CHASE - confirm data type is not list.@

    :param sources: @CHASE - This parameter is not used and is just passed back,
    can be removed from the call?@

    :type sources: numpy.ndarray @CHASE - Is this a list/numpy.ndarray?@
    """

    ras, decs, responses = bg_mask_annulus(band, ra0, dec0, annulus, ras,
                                           decs, responses)

    return bg_mask_sources(band, ra0, dec0, ras, decs, responses, sources)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def cheese_bg_area(band, ra0, dec0, annulus, sources, nsamples=10e5, ntests=10):
    """
    @CHASE - Please provide descriptor.@

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: Right ascension, in degrees, of the center of the aperture.
    @CHASE - confirm please.@

    :type ra0: float

    :param dec0: Declination, in degrees, of the center of the aperture.
    @CHASE - confirm please.@

    :type dec0: float

    :param annulus: The size of the inner and outer background annuli,
    in degrees. @CHASE - please confirm units.@

    :type annulus: list

    :param sources: Set of sources to check whether they are in the background.

    :type sources: numpy.ndarray @CHASE - Is this a list or numpy.ndarray?@

    :param nsamples: Number of trials to run per test.

    :type nsamples: int

    :param ntests: Number of tests to run.

    :type ntests: int

    :returns: numpy.ndarray - The size of the annuli for each trial.
    """

    # This is just a really naive Monte Carlo.
    ratios = np.zeros(ntests)

    for i in range(ntests):
        ann_events = bg_mask_annulus(band, ra0, dec0, annulus,
                                     np.random.uniform(ra0-annulus[1],
                                                       ra0+annulus[1],
                                                       int(nsamples)),
                                     np.random.uniform(dec0-annulus[1],
                                                       dec0+annulus[1],
                                                       int(nsamples)),
                                     np.ones(nsamples))
        mask_events = bg_mask_sources(band, ra0, dec0, ann_events[0],
                                      ann_events[1], ann_events[2], sources)

        try:
            ratios[i] = float(mask_events[2].sum())/float(ann_events[2].sum())
        except ZeroDivisionError:
            ratios[i] = 0.

    return (mc.area(annulus[1])-mc.area(annulus[0]))*ratios.mean()
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def cheese_bg(band, ra0, dec0, radius, annulus, ras, decs, responses,
              maskdepth=20., maskradius=1.5, eff_area=False, sources=False):
    """
    Returns an estimate of the number of counts (not count rate) within the
    aperture based upon a masked background annulus.

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: Right ascension, in degrees, of the center of the aperture.
    @CHASE - confirm please.@

    :type ra0: float

    :param dec0: Declination, in degrees, of the center of the aperture.
    @CHASE - confirm please.@

    :type dec0: float

    :param radius: The photometric aperture, in degrees. @CHASE - confirm
    units.@

    :type radius: float

    :param annulus: The size of the inner and outer background annuli,
    in degrees. @CHASE - please confirm units.@

    :type annulus: list

    :param ras: Right ascension of sources to check whether they lie in the
    annulus, in degrees.

    :type ras: numpy.ndarray @CHASE - confirm data type is not list.@

    :param decs: Declination of sources to check whether they lie in the
    annulus, in degrees.

    :type decs: numpy.ndarray @CHASE - confirm data type is not list.@

    :param responses: Set of response values for the sources that will be
    checked whether they lie in the annulus.

    :type responses: numpy.ndarray @CHASE - confirm data type is not list.@

    :param maskdepth: @CHASE - please provide description.@

    :type maskdepth: float

    :param maskradius: @CHASE - please provide description.@

    :type maskradius: float

    :param eff_area: @CHASE - please provide description.@

    :type eff_area: bool

    :param sources: @CHASE - please provide description.@

    :type sources: bool @CHASE - This looks like a mixed-type, the default
    should probably be None and data type list or numpy.ndarray?@

    :returns: float - The number of counts (excluding background) within the
    photometric apreture. @CHASE - Confirm description, this does not count
    background?@
    """

    # [Future]: This recomputes eff_area every pass at huge computational cost.
    if not sources:
        sources = bg_sources(band, ra0, dec0, annulus[1], maskdepth=maskdepth)

    bg_counts = bg_mask(band, ra0, dec0, annulus, ras, decs, responses,
                        sources)[2].sum()
    if not eff_area:
        eff_area = cheese_bg_area(band, ra0, dec0, annulus, sources)

    return mc.area(radius)*bg_counts/eff_area if eff_area else 0.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def reduce_lcurve(bin_ix, region_ix, data, function, dtype='float64'):
    """
    Produces light curve columns by iteratively applying 'function' to 'data'
    within 'region_ix' over 'bin_ix'.

    :param bin_ix: @CHASE - please describe parameter.@

    :type bin_ix: list @CHASE - confirm data type please.@

    :param region_ix: @CHASE - please describe parameter.@

    :type region_ix: list @CHASE - confirm data type please.@

    :param data: The data to apply the function on.

    :type data: numpy.ndarray @CHASE - list or numpy.ndarray?@

    :param function: The function to apply to the data.

    :type function: function

    :param dtype: The data type.

    :type dtype: str

    :returns: numpy.ndarray -- The light curve columns.
    """

    bin_num = np.unique(bin_ix)
    output = np.empty(len(bin_num))

    for i, b in enumerate(bin_num):
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
def maskwarning(band, bin_ix, events, verbose=0, mapkey='H'):
    """
    Test if any given events are near a masked detector region.

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param bin_ix: @CHASE - please provide description.@

    :type bin_ix: numpy.ndarray @CHASE - is this list or numpy.ndarray?@

    :param events: Set of photon events to check if they are near a masked
    detector region.

    :type events: dict

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param mapkey: @CHASE - please provide description.@

    :type mapkey: str

    :returns: bool -- Returns True/False whether a given set of events are near
    a masked detector region.
    """

    maps = {'H':cal.mask, 'E':cal.flat}

    img, _ = maps[mapkey](band, buffer=True)

    ix = np.where(
        img[np.array(events['photons']['col'][bin_ix], dtype='int16'),
            np.array(events['photons']['row'][bin_ix], dtype='int16')] == 0)

    return True if len(ix[0]) else False
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def lowresponsewarning(bin_ix, events, verbose=0, ratio=0.7):
    """
    @CHASE - Please provide descriptor.@

    :param bin_ix: @CHASE - please provide description.@

    :type bin_ix: numpy.ndarray @CHASE - is this list or numpy.ndarray?@

    :param events: Set of photon events to check if there is a low response.

    :type events: dict

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param ratio: The value that defines a low response, between 0 and 1.
    @CHASE - please refine, this is a response 'ratio', and if so, what is
    the ratio?@

    :type ratio: float

    :returns: bool -- Returns True/False whether a given set of events have
    a low response. @CHASE - please refine description.@
    """

    ix = np.where(events['photons']['response'][bin_ix] < 0.7)

    return True if len(ix[0]) else False
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def exptimewarning(bin_ix, events, verbose=0, ratio=0.5):
    """
    @CHASE - Please provide descriptor.@

    :param bin_ix: @CHASE - please provide description.@

    :type bin_ix: numpy.ndarray @CHASE - is this list or numpy.ndarray?@

    :param events: Set of photon events to check if there is an effective
    exposure time warning.

    :type events: dict

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param ratio: The value that defines an exposure time warning, between 0
    and 1.
    @CHASE - please refine, what is this a ratio of?@

    :type ratio: float

    :returns: bool -- Returns True/False whether a given set of events have
    a low effective exposure time. @CHASE - please refine description.@
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

    :param bin_ix: @CHASE - please provide description.@

    :type bin_ix: numpy.ndarray @CHASE - is this list or numpy.ndarray?@

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
    @CHASE - Please provide descriptor.@

    :param bin_ix: @CHASE - please provide description.@

    :type bin_ix: numpy.ndarray @CHASE - is this list or numpy.ndarray?@

    :param events: Set of photon events to check if they are near the detector
    edge.

    :type events: dict

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param valid_detrad: The radius, in degrees, beyond which an edge warning is
    raised.

    :type detrad: float

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
    1 - 'hotspot' - events in pixels contiguous to a hotspot masked region
    2 - 'mask edge' - events in pixels contiguous to the detector edge
    4 - 'exptime' - bin contains < 50% exposure time coverage
    8 - 'respose' - events weighted with response < 0.7
    16 - 'nonlinearity' - local countrate exceeds 10% response dropoff
    32 - 'detector edge' - events outside of 0.5 degrees of detector center

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param bin_ix: @CHASE - please provide description.@

    :type bin_ix: numpy.ndarray @CHASE - is this list or numpy.ndarray?@

    :param events: Set of photon events to check for warning flags.

    :type events: dict

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :returns: numpy.ndarray -- The array of flags for each photon event.
    """

    bin_num = np.unique(bin_ix)
    flags = np.zeros(len(bin_num))

    for i, b in enumerate(bin_num):
        try:
            ix = bin_ix[np.where(bin_ix == bin)]

            if maskwarning(band, ix, events, mapkey='H', verbose=verbose):
                flags[i] += 1
            if maskwarning(band, ix, events, mapkey='E', verbose=verbose):
                flags[i] += 2
            if exptimewarning(i, events, verbose=verbose):
                flags[i] += 4
            if lowresponsewarning(ix, events, verbose=verbose):
                flags[i] += 8
            if nonlinearitywarning(band, i, events, verbose=verbose):
                flags[i] += 16
            if detedgewarning(ix, events, verbose=verbose):
                flags[i] += 32

        except IndexError:
            return np.array([np.nan])
        except:
            raise

    return np.array(flags)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def quickmag(band, ra0, dec0, tranges, radius, annulus=None, data={},
             stepsz=None, verbose=0, detsize=1.25, coadd=False,
             photonfile=None):
    """
    @CHASE - Please provide descriptor.@

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: Right ascension, in degrees, of the center of the field-of-view.
    @CHASE - confirm please.@

    :type ra0: float

    :param dec0: Declination, in degrees, of the center of the field-of-view.
    @CHASE - confirm please.@

    :type dec0: float

    :param tranges: Set of time ranges to query within. @CHASE - GALEX time?@

    :type tranges: list

    :param radius: The photometric aperture, in degrees. @CHASE - Confirm
    this is in degrees.@

    :type radius: float

    :param annulus: Radius of the inner and outer annuli to define the
    background with, in degrees. @CHASE - confirm units.@

    :type annulus: list

    :param data: Set of photon events to use. @CHASE - Is this dict updated
    within this method or sub-methods? Otherwise, why is it provided on input?@

    :type data: dict

    :param stepsz: The size of the time bins to use, in seconds.

    :type stepsz: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :param coadd: Set to True if calculating a total flux instead of flux
    from each time bin.

    :type coadd: bool

    :param photonfile: Name of photon event CSV file to use. @CHASE - This is
    not used in this method, it can be removed?@

    :type photonfile: str

    :returns: dict -- The light curve, including input parameters.
    """

    if verbose:
        mc.print_inline("Retrieving all of the target events.")
    searchradius = radius if annulus is None else annulus[1]
    data = pullphotons(band, ra0, dec0, tranges, searchradius, verbose=verbose)

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
            dbt.compute_exptime(band, tranges if coadd else zip(lcurve['t0'],
                                                                lcurve['t1']),
                                verbose=verbose, coadd=coadd, detsize=detsize,
                                skypos=[ra0, dec0]))
    except IndexError:
        if np.isnan(data['t']):
            print "No valid data available in {t}".format(t=tranges)
        lcurve['t0'] = np.array([np.nan])
        lcurve['t1'] = np.array([np.nan])
        lcurve['exptime'] = np.array([0])

    angSep = mc.angularSeparation(ra0, dec0, data['ra'], data['dec'])

    aper_ix = np.where(angSep <= radius)
    lcurve['t0_data'] = reduce_lcurve(bin_ix, aper_ix, data['t'], np.min)
    lcurve['t1_data'] = reduce_lcurve(bin_ix, aper_ix, data['t'], np.max)
    lcurve['t_mean'] = reduce_lcurve(bin_ix, aper_ix, data['t'], np.mean)
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

    lcurve['mcat_bg'] = lcurve['exptime']*np.array(
        [dbt.mcat_skybg(band, [ra0, dec0], radius, trange=tr, verbose=verbose)
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
def getcurve(band, ra0, dec0, radius, annulus=None, stepsz=None, lcurve={},
             trange=None, tranges=None, verbose=0, coadd=False, minexp=1.,
             maxgap=1., photonfile=None, detsize=1.1):
    """
    @CHASE - Please provide descriptor. Also, although maybe not worth
    tacking now, this should probably be called 'get_curve' since the sibling
    method is called 'write_curve'.@

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: Right ascension, in degrees, of the center of the field-of-view.
    @CHASE - confirm please.@

    :type ra0: float

    :param dec0: Declination, in degrees, of the center of the field-of-view.
    @CHASE - confirm please.@

    :type dec0: float

    :param radius: The photometric aperture, in degrees. @CHASE - Confirm
    this is in degrees.@

    :type radius: float

    :param annulus: Radius of the inner and outer annuli to define the
    background with, in degrees. @CHASE - confirm units.@

    :type annulus: list

    :param stepsz: The size of the time bins to use, in seconds.

    :type stepsz: float

    :param lcurve: @CHASE - why is this passed as an argument?@

    :type lcurve: dict

    :param trange: Minimum and maximum time range to make a light curve.
    @CHASE - assume this is in GALEX time?@

    :type trange: list

    :param tranges: Set of time ranges to query within. @CHASE - GALEX time?@

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

    :param photonfile: Name of photon event CSV file to use. @CHASE - This is
    not used in this method, it can be removed?@

    :type photonfile: str

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :returns: dict -- The light curve, including input parameters.
    """

    skyrange = [np.array(annulus).max().tolist() if annulus else radius,
                np.array(annulus).max().tolist() if annulus else radius,]
    if verbose:
        mc.print_inline("Getting exposure ranges.")
    if tranges is None:
        tranges = dbt.fGetTimeRanges(band, [ra0, dec0], trange=trange,
                                     maxgap=maxgap, minexp=minexp,
                                     verbose=verbose, detsize=detsize)
    if not np.array(tranges).shape[1]:
        print "No exposure time at this location: [{ra},{dec}]".format(
            ra=ra0, dec=dec0)
        return None
    lcurve = quickmag(band, ra0, dec0, tranges, radius, annulus=annulus,
                      stepsz=stepsz, verbose=verbose, coadd=coadd)

    if verbose:
        mc.print_inline("Done.")
        mc.print_inline("")

    return lcurve
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def write_curve(band, ra0, dec0, radius, csvfile=None, annulus=None,
                stepsz=None, trange=None, tranges=None, verbose=0, coadd=False,
                iocode='wb', detsize=1.1, overwrite=False, minexp=1., maxgap=1.,
                photonfile=None, minimal_output=False):
    """
    @CHASE - Please provide descriptor.@

    :param band: The band being used, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: Right ascension, in degrees, of the center of the field-of-view.
    @CHASE - confirm please.@

    :type ra0: float

    :param dec0: Declination, in degrees, of the center of the field-of-view.
    @CHASE - confirm please.@

    :type dec0: float

    :param radius: The photometric aperture, in degrees. @CHASE - Confirm
    this is in degrees.@

    :type radius: float

    :param csvfile: Name of the photon event CSV file to use.

    :type csvfile: str

    :param annulus: Radius of the inner and outer annuli to define the
    background with, in degrees. @CHASE - confirm units.@

    :type annulus: list

    :param stepsz: The size of the time bins to use, in seconds.

    :type stepsz: float

    :param trange: Minimum and maximum time range to make a light curve.
    @CHASE - assume this is in GALEX time?@

    :type trange: list

    :param tranges: Set of time ranges to query within. @CHASE - GALEX time?@

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

    :param photonfile: Name of photon event CSV file to use. @CHASE - This is
    not used in this method or in 'getcurve', it can be removed?@

    :type photonfile: str

    :param minimal_output: If True, produce an output file with a minimum
    number of columns.

    :type minimal_output: bool

    :returns: dict -- The light curve, including input parameters.
    """

    data = getcurve(band, ra0, dec0, radius, annulus=annulus, stepsz=stepsz,
                    trange=trange, tranges=tranges, verbose=verbose,
                    coadd=coadd, minexp=minexp, maxgap=maxgap,
                    photonfile=photonfile, detsize=detsize)
    if not data:
        return None

    exclude_keys = ['photons', 'params']
    minimal_columns = ['t0', 't1', 'exptime', 'flux', 'flux_err', 'flags']

    if annulus:
        minimal_columns += ['flux_bgsub', 'flux_bgsub_err']

    if csvfile:
        if verbose:
            mc.print_inline('Building output data frame.')

        frame, columns = {}, []

        for k in minimal_columns if minimal_output else data.keys():
            if k in exclude_keys:
                continue
            frame[k] = data[k]
            columns += [k]

        try:
            output = pd.DataFrame(frame)
        except:
            raise
            if verbose > 1:
                print 'Unable to build dataframe.'
        try:
            output.to_csv(csvfile, index=False, mode=iocode, columns=columns)
        except:
            print 'Unable to write to: '+str(csvfile)

    else:
        if verbose > 2:
            print "No CSV file requested."

        if verbose or (not verbose and not csvfile):
            print "AB Magnitudes:               "
            print data['mag']

    return data
# ------------------------------------------------------------------------------
