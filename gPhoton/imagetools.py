"""
.. module:: imagetools

   :synopsis: @CHASE - please provide summary of this module.@

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

import gQuery
import numpy as np
import MCUtils as mc
from astropy import wcs as pywcs
from astropy.io import fits as pyfits
import scipy.misc
import scipy.special # erfc
import scipy.ndimage
import dbasetools as dbt
import galextools as gxt
import curvetools as ct
from gQuery import tscale
from gPhoton import __version__

# ------------------------------------------------------------------------------
def define_wcs(skypos, skyrange, width=False, height=False, verbose=0,
               pixsz=0.000416666666666667):
    """
    Define the world coordinate system (WCS).

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list @CHASE - list or numpy.ndarray?@

    :param skyrange: @CHASE - please describe.@

    :type skyrange: list @CHASE - list or numpy.ndarray?@

    :param width: @CHASE - This is not used, can be removed?@

    :type width: bool

    :param height: @CHASE - This is not used, can be removed?@

    :type height: bool

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param pixsz: Size of a GALEX pixel, in degrees. @CHASE - confirm unit@

    :type pixsz: float

    :returns: astropy.wcs WCS Object -- The WCS information.
    """

    if verbose:
        mc.print_inline('Defining World Coordinate System (WCS).')

    # NAXIS = 2
    wcs = pywcs.WCS(naxis=2)

    imsz = gxt.deg2pix(skypos, skyrange)

    wcs.wcs.cdelt = np.array([-pixsz, pixsz])
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    wcs.wcs.crpix = [(imsz[1]/2.)+0.5, (imsz[0]/2.)+0.5]
    wcs.wcs.crval = skypos

    return wcs
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def movie_tbl(band, tranges, verbose=0, framesz=0, retries=20):
    """
    Initialize a FITS table to contain movie frame information.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param tranges: Set of time ranges to retrieve the photon events.
    @CHASE - in GALEX time?@

    :type tranges: list @CHASE - list or numpy.ndarray?@

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param framesz: The time bin size to use per frame, in seconds. @CHASE -
    please confirm description.@

    :type framesz: int @CHASE - needs to be int or can be float?@

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int

    :returns: astropy.fits.BinTableHDU object -- The set of frames as an HDU
    object. @CHASE - check return object type.@
    """

    if verbose:
        mc.print_inline('Populating exposure time table.')

    tstarts, tstops, exptimes = [], [], []

    for trange in tranges:
        stepsz = framesz if framesz else trange[1]-trange[0]
        steps = np.ceil((trange[1]-trange[0])/stepsz)
        for i, t0 in enumerate(np.arange(trange[0], trange[1], stepsz)):
            t1 = trange[1] if i == steps else t0+stepsz
            tstarts.append(t0)
            tstops.append(t1)
            exptimes.append(dbt.compute_exptime(band, [t0, t1],
                                                verbose=verbose,
                                                retries=retries))
    col1 = pyfits.Column(name='tstart', format='E', array=np.array(tstarts))
    col2 = pyfits.Column(name='tstop', format='E', array=np.array(tstops))
    col3 = pyfits.Column(name='exptime', format='E', array=np.array(exptimes))
    cols = pyfits.ColDefs([col1, col2, col3])
    tbl = pyfits.BinTableHDU.from_columns(cols)

    return tbl
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def fits_header(band, skypos, tranges, skyrange, width=False, height=False,
                verbose=0, hdu=False, retries=20):
    """
    Populate a FITS header.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list @CHASE - list or numpy.ndarray?@

    :param tranges: Set of time ranges to retrieve the photon events.
    @CHASE - in GALEX time?@

    :type tranges: list @CHASE - list or numpy.ndarray?@

    :param skyrange: @CHASE - please describe.@

    :type skyrange: list @CHASE - list or numpy.ndarray?@

    :param width: @CHASE - This is not used in define_wcs, can be removed?@

    :type width: bool

    :param height: @CHASE - This is not used in define_wcs, can be removed?@

    :type height: bool

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param hdu: An existing HDU to modify. @CHASE - should default should be
    None instead of a bool?@

    :type hdu: bool

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int

    :returns: astropy.fits.BinTableHDU object -- The modified HDU header.
    @CHASE - check return object type and refine description if needed.@
    """

    if verbose:
        mc.print_inline('Populating FITS header.')

    hdu = hdu if hdu else pyfits.PrimaryHDU()

    wcs = define_wcs(skypos, skyrange, width=width, height=height)

    hdu.header['CDELT1'], hdu.header['CDELT2'] = wcs.wcs.cdelt
    hdu.header['CTYPE1'], hdu.header['CTYPE2'] = wcs.wcs.ctype
    hdu.header['CRPIX1'], hdu.header['CRPIX2'] = wcs.wcs.crpix
    hdu.header['CRVAL1'], hdu.header['CRVAL2'] = wcs.wcs.crval
    hdu.header['EQUINOX'], hdu.header['EPOCH'] = 2000., 2000.
    hdu.header['BAND'] = 1 if band == 'NUV' else 2
    hdu.header['VERSION'] = 'v{v}'.format(v=__version__)

    return hdu
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def makemap(band, skypos, trange, skyrange, response=False, verbose=0,
            detsize=1.1):
    """
    @CHASE - please provide description of this method.@
    @CHASE - If remove 'width' and 'height' from define_wcs, make sure to remove
    in the call here.@

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list @CHASE - list or numpy.ndarray?@

    :param trange: Minimum and maximum time to use, in GALEX time.

    :type trange: list @CHASE - list or numpy.ndarray?@

    :param skyrange: @CHASE - please describe.@

    :type skyrange: list @CHASE - list or numpy.ndarray?@

    :param response: @CHASE - please explain what this option does.@

    :type response: bool

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :returns: numpy.ndarray - The bi-dimensional histogram of x and y.
    """

    imsz = gxt.deg2pix(skypos, skyrange)

    photons = np.array(gQuery.getArray(
        gQuery.skyrect(band, skypos[0], skypos[1], trange[0], trange[1],
                       skyrange[0], skyrange[1]), verbose=verbose),
                       dtype='float64')
    try:
        events = {'t':photons[:, 0]/tscale, 'ra':photons[:, 1],
                  'dec':photons[:, 2], 'xi':photons[:, 3], 'eta':photons[:, 4],
                  'x':photons[:, 5], 'y':photons[:, 6]}
    except IndexError:
        if verbose > 2:
            print ('No events found at {s} +/- {r} in {t}.'.format(
                s=skypos, r=skyrange, t=trange))
        return np.zeros(imsz)

    # Trim the data on detsize
    col, row = ct.xieta2colrow(events['xi'], events['eta'], band)
    ix = np.where((1.25/800.)*mc.distance(col, row, 400, 400) <= detsize)
    n = len(ix[0])
    m = len(col)

    if n == 0:
        return np.zeros(imsz)

    for k in events.keys():
        events[k] = events[k][ix]

    events = ct.hashresponse(band, events)
    wcs = define_wcs(skypos, skyrange, width=False, height=False)
    coo = zip(events['ra'], events['dec'])
    foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo, 1), 1)
    weights = 1./events['response'] if response else None
    H, xedges, yedges = np.histogram2d(foc[:, 1]-0.5, foc[:, 0]-0.5, bins=imsz,
                                       range=([[0, imsz[0]], [0, imsz[1]]]),
                                       weights=weights)

    return H
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def integrate_map(band, skypos, tranges, skyrange, width=False, height=False,
                  verbose=0, memlight=False, hdu=False, retries=20,
                  response=False, detsize=1.1):
    """
    Integrate an image over some number of time ranges. Use a reduced
	memory optimization (at the expense of more web queries) if requested.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list @CHASE - list or numpy.ndarray?@

    :param tranges: Set of time ranges to retrieve the photon events.
    @CHASE - in GALEX time?@

    :type tranges: list @CHASE - list or numpy.ndarray?@

    :param skyrange: @CHASE - please describe.@

    :type skyrange: list @CHASE - list or numpy.ndarray?@

    :param width: @CHASE - This is not used in called funct., can be removed?@

    :type width: bool

    :param height: @CHASE - This is not used in called funct., can be removed?@

    :type height: bool

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param memlight: Reduce memory usage by breaking query into smaller
    segments. @CHASE - refine if needed, is this the size in seconds?@

    :type memlight: bool @CHASE - Should be None instead of bool?@

    :param hdu: An existing HDU to modify. @CHASE - should default should be
    None instead of a bool?@

    :type hdu: bool

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int

    :param response: @CHASE - please explain what this option does.@

    :type response: bool

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :returns: numpy.ndarray - The integrated image.
    """

    imsz = gxt.deg2pix(skypos, skyrange)
    img = np.zeros(imsz)

    for trange in tranges:
        # If memlight is requested, break the integration into
        # smaller chunks.
        # [Future]: memlight gives slightly wrong answers right now
        # This is probably due to a quirk of SQL, per issue #140.
        # Deprecating memlight until this can be resolved.
        step = memlight if memlight else trange[1]-trange[0]

        for i in np.arange(trange[0], trange[1], step):
            t0, t1 = i, i+step
            if verbose:
                mc.print_inline('Coadding '+str(t0)+' to '+str(t1))
            img += makemap(band, skypos, [t0, t1], skyrange, response=response,
                           verbose=verbose, detsize=detsize)

        if response: # This is an intensity map.
            img /= dbt.compute_exptime(band, trange, skypos=skypos,
                                       verbose=verbose)

    return img
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def write_jpeg(filename, band, skypos, tranges, skyrange, width=False,
               height=False, stepsz=1., overwrite=False, verbose=0, retries=20):
    """
    Write a 'preview' jpeg image from a count map.

    :param filename: Name of output jpeg to make.

    :type filename: str

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list @CHASE - list or numpy.ndarray?@

    :param tranges: Set of time ranges to retrieve the photon events.
    @CHASE - in GALEX time?@

    :type tranges: list @CHASE - list or numpy.ndarray?@

    :param skyrange: @CHASE - please describe.@

    :type skyrange: list @CHASE - list or numpy.ndarray?@

    :param width: @CHASE - This is not used in called funct., can be removed?@

    :type width: bool

    :param height: @CHASE - This is not used in called funct., can be removed?@

    :type height: bool

    :param stepsz: Time bin size to use, in seconds.

    :type stepsz: float

    :param overwrite: Overwrite existing output files?

    :type overwrite: bool

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int
    """

    scipy.misc.imsave(filename, integrate_map(band, skypos, tranges, skyrange,
                                              width=width, height=height,
                                              verbose=verbose, retries=retries))

    # @CHASE - No need for this return?@
    return
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def movie(band, skypos, tranges, skyrange, framesz=0, width=False, height=False,
          verbose=0, memlight=False, coadd=False, response=False, hdu=False,
          retries=20, detsize=1.1):
    """
    Generate a movie (mov) file.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list @CHASE - list or numpy.ndarray?@

    :param tranges: Set of time ranges to retrieve the photon events.
    @CHASE - in GALEX time?@

    :type tranges: list @CHASE - list or numpy.ndarray?@

    :param skyrange: @CHASE - please describe.@

    :type skyrange: list @CHASE - list or numpy.ndarray?@

    :param framesz: The time bin size to use per frame, in seconds. @CHASE -
    please confirm description.@

    :type framesz: int @CHASE - needs to be int or can be float?@

    :param width: @CHASE - This is not used in called funct., can be removed?@

    :type width: bool

    :param height: @CHASE - This is not used in called funct., can be removed?@

    :type height: bool

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param memlight: Reduce memory usage by breaking query into smaller
    segments. @CHASE - refine if needed, is this the size in seconds?@

    :type memlight: bool @CHASE - Should be None instead of bool?@

    :param coadd: Create a coadd movie. @CHASE - can you provide better
    description.@

    :type coadd: bool

    :param response: @CHASE - please explain what this option does.@

    :type response: bool

    :param hdu: An existing HDU to modify. @CHASE - should default should be
    None instead of a bool?@

    :type hdu: bool

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :returns: numpy.ndarray -- The movie file.
    """

    # Not defining stepsz creates a single full depth image.
    if verbose:
        print tranges

    if coadd or (len(tranges) == 1 and not framesz) or (not len(tranges)):
        if verbose > 2:
            print 'Coadding across '+str(tranges)

        mv = integrate_map(band, skypos, tranges, skyrange, width=width,
                           height=height, verbose=verbose, memlight=memlight,
                           hdu=hdu, retries=retries, response=response,
                           detsize=detsize)
    else:
        for trange in tranges:
            stepsz = framesz if framesz else trange[1]-trange[0]
            steps = np.ceil((trange[1]-trange[0])/stepsz)
            for i, t0 in enumerate(np.arange(trange[0], trange[1], stepsz)):
                if verbose > 1:
                    mc.print_inline('Movie frame '+str(i+1)+' of '+
                                    str(int(steps)))
                t1 = trange[1] if i == steps else t0+stepsz

                img = integrate_map(band, skypos, [[t0, t1]], skyrange,
                                    width=width, height=height, verbose=verbose,
                                    memlight=memlight, hdu=hdu, retries=retries,
                                    response=response, detsize=detsize)
                if img.min() == 0 and img.max() == 0:
                    if verbose > 1:
                        print 'No data in frame {i}. Skipping...'.format(i=i)
                    continue
                try:
                    mv.append(img)
                except:
                    mv = [img]

    return np.array(mv)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def create_image(band, skypos, tranges, skyrange, framesz=0, width=False,
                 height=False, verbose=0, memlight=False, coadd=False,
                 response=False, hdu=False, retries=20, detsize=1.1):
    """
    @CHASE - please provide description of this method.@

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list @CHASE - list or numpy.ndarray?@

    :param tranges: Set of time ranges to retrieve the photon events.
    @CHASE - in GALEX time?@

    :type tranges: list @CHASE - list or numpy.ndarray?@

    :param skyrange: @CHASE - please describe.@

    :type skyrange: list @CHASE - list or numpy.ndarray?@

    :param framesz: The time bin size to use per frame, in seconds. @CHASE -
    please confirm description.@

    :type framesz: int @CHASE - needs to be int or can be float?@

    :param width: @CHASE - This is not used in called funct., can be removed?@

    :type width: bool

    :param height: @CHASE - This is not used in called funct., can be removed?@

    :type height: bool

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param memlight: Reduce memory usage by breaking query into smaller
    segments. @CHASE - refine if needed, is this the size in seconds?@

    :type memlight: bool @CHASE - Should be None instead of bool?@

    :param coadd: Create a coadd movie. @CHASE - can you provide better
    description.@

    :type coadd: bool

    :param response: @CHASE - please explain what this option does.@

    :type response: bool

    :param hdu: An existing HDU to modify. @CHASE - should default should be
    None instead of a bool?@

    :type hdu: bool

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :returns: numpy.ndarray -- The image file.
    """

    img = movie(band, skypos, tranges, skyrange, framesz=framesz,
                width=width, height=height, verbose=verbose, memlight=memlight,
                coadd=coadd, response=response, hdu=hdu, retries=retries,
                detsize=detsize)

    return np.array(img)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def write_images(band, skypos, tranges, skyrange, write_cnt=False,
                 write_int=False, write_rr=False, framesz=0, width=False,
                 height=False, verbose=0, memlight=False, coadd=False,
                 overwrite=False, retries=20, write_cnt_coadd=False,
                 write_int_coadd=False, detsize=1.1):
    """
    Generate a write various maps to files.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list @CHASE - list or numpy.ndarray?@

    :param tranges: Set of time ranges to retrieve the photon events.
    @CHASE - in GALEX time?@

    :type tranges: list @CHASE - list or numpy.ndarray?@

    :param skyrange: @CHASE - please describe.@

    :type skyrange: list @CHASE - list or numpy.ndarray?@

    :param write_cnt: Make count image?

    :type write_cnt: bool

    :param write_int: Make intensity image?

    :type write_int: bool

    :param write_rr: Make relative response image?
    @CHASE - Not supported/used any more it looks like, since it was commented
    out.  I removed for clarity, should remove from options (for now?)

    :type write_rr: bool

    :param framesz: The time bin size to use per frame, in seconds. @CHASE -
    please confirm description.@

    :type framesz: int @CHASE - needs to be int or can be float?@

    :param width: @CHASE - This is not used in called funct., can be removed?@

    :type width: bool

    :param height: @CHASE - This is not used in called funct., can be removed?@

    :type height: bool

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param memlight: Reduce memory usage by breaking query into smaller
    segments. @CHASE - refine if needed, is this the size in seconds?@

    :type memlight: bool @CHASE - Should be None instead of bool?@

    :param coadd: Create a coadd movie. @CHASE - can you provide better
    description.@

    :type coadd: bool

    :param overwrite: Overwrite existing output files?

    :type overwrite: bool

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int

    :param write_cnt_coadd: Make count coadd image?

    :type write_cnt_coadd: bool

    :param write_int_coadd: Make intensity coadd image?

    :type write_int_coadd: bool

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float
    """

    # No files were requested, so don't bother doing anything.
    imtypes = {'cnt':write_cnt, 'int':write_int, 'int_coadd':write_int_coadd,
               'cnt_coadd':write_cnt_coadd}

    for i in imtypes.keys():
        if not imtypes[i]:
            continue

        img = create_image(band, skypos, tranges, skyrange, framesz=framesz,
                           width=width, height=height, verbose=verbose,
                           memlight=memlight, retries=retries, detsize=detsize,
                           coadd=(
                               True if (coadd or i in ['cnt_coadd',
                                                       'int_coadd']) else
                               False),
                           response=(
                               True if i in ['int', 'int_coadd'] else False))

        # Add a conditional so that this is only created for multi-frame images
        tbl = movie_tbl(band, tranges, framesz=framesz, verbose=verbose,
                        retries=retries) if i in ['int', 'int_coadd'] else False

        hdu = pyfits.PrimaryHDU(img)
        hdu = fits_header(band, skypos, tranges, skyrange, width=width,
                          height=height, verbose=verbose, hdu=hdu,
                          retries=retries)

        hdulist = pyfits.HDUList([hdu, tbl]) if tbl else pyfits.HDUList([hdu])

        if verbose:
            print 'Writing image to {o}'.format(o=imtypes[i])

        hdulist.writeto(imtypes[i], clobber=overwrite)

    # @CHASE - this return is not needed?@
    return
# ------------------------------------------------------------------------------
