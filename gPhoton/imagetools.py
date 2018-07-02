"""
.. module:: imagetools
   :synopsis: Tools for the cration of count and intensity images and movies.
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
from astropy import wcs as pywcs
from astropy.io import fits as pyfits
from builtins import str
from builtins import zip
import numpy as np
import scipy.misc
import scipy.ndimage
import scipy.special # erfc
# gPhoton imports.
import gPhoton.curvetools as ct
import gPhoton.dbasetools as dbt
import gPhoton.galextools as gxt
from gPhoton import __version__
import gPhoton.MCUtils as mc
import gPhoton.gQuery as gQuery
from gPhoton.gQuery import tscale

# ------------------------------------------------------------------------------
def define_wcs(skypos, skyrange, verbose=0, pixsz=0.000416666666666667):
    """
    Define the world coordinate system (WCS).

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param skyrange: Extent of the region of interest, in degrees.

    :type skyrange: list

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param pixsz: Size of a GALEX pixel, in degrees.

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
def movie_tbl(band, tranges, verbose=0, framesz=0., retries=100):
    """
    Initialize a FITS table to contain movie frame information.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param tranges: Set of time ranges to retrieve the photon events, in
        GALEX time seconds.

    :type tranges: list

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param framesz: The time bin size (depth)to use per frame, in seconds.

    :type framesz: float

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int

    :returns: astropy.fits.BinTableHDU object -- The set of frames as an HDU
        object.
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
                                                verbose=verbose))
    col1 = pyfits.Column(name='tstart', format='E', array=np.array(tstarts))
    col2 = pyfits.Column(name='tstop', format='E', array=np.array(tstops))
    col3 = pyfits.Column(name='exptime', format='E', array=np.array(exptimes))
    cols = pyfits.ColDefs([col1, col2, col3])
    tbl = pyfits.BinTableHDU.from_columns(cols)

    return tbl
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def fits_header(band, skypos, tranges, skyrange, verbose=0, hdu=None,
                retries=100):
    """
    Populate a FITS header.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param tranges: Set of time ranges to retrieve the photon events, in
        GALEX time seconds.

    :type tranges: list

    :param skyrange: Extent of the region of interest, in degrees.

    :type skyrange: list

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param hdu: An existing HDU to modify.

    :type hdu: bool

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int

    :returns: astropy.fits.BinTableHDU object -- The modified HDU header.
    """

    if verbose:
        mc.print_inline('Populating FITS header.')

    hdu = hdu if hdu else pyfits.PrimaryHDU()

    wcs = define_wcs(skypos, skyrange)

    hdu.header['CDELT1'], hdu.header['CDELT2'] = wcs.wcs.cdelt
    hdu.header['CTYPE1'], hdu.header['CTYPE2'] = wcs.wcs.ctype
    hdu.header['CRPIX1'], hdu.header['CRPIX2'] = wcs.wcs.crpix
    hdu.header['CRVAL1'], hdu.header['CRVAL2'] = wcs.wcs.crval
    hdu.header['EQUINOX'], hdu.header['EPOCH'] = 2000., 2000.
    hdu.header['BAND'] = 1 if band == 'NUV' else 2
    hdu.header['VERSION'] = 'v{v}'.format(v=__version__)
    hdu.header['EXPSTART'] = np.array(tranges).min()
    hdu.header['EXPEND'] = np.array(tranges).max()
    hdu.header['EXPTIME'] = sum(t1-t0 for (t0,t1) in tranges)

    return hdu
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def makemap(band, skypos, trange, skyrange, response=False, verbose=0,
            detsize=1.1):
    """
    Generate a single image frame.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param trange: Minimum and maximum time to use, in GALEX time seconds.

    :type trange: list

    :param skyrange: RA and Dec extent of the region of interest in degrees.

    :type skyrange: list

    :param response: Apply the response correction.

    :type response: bool

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :returns: numpy.ndarray - The bi-dimensional histogram of ra and dec.
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
            print('No events found at {s} +/- {r} in {t}.'.format(
                s=skypos, r=skyrange, t=trange))
        return np.zeros(np.array(imsz, dtype='int32'))

    # Trim the data on detsize
    col, row = ct.xieta2colrow(events['xi'], events['eta'], band)
    ix = np.where(gxt.aper2deg(4)*mc.distance(col, row, 400, 400) <= detsize)
    n = len(ix[0])
    m = len(col)

    if n == 0:
        return np.zeros(np.int(imsz))

    for k in list(events.keys()):
        events[k] = events[k][ix]

    events = ct.hashresponse(band, events)
    wcs = define_wcs(skypos, skyrange)
    coo = list(zip(events['ra'], events['dec']))
    foc = wcs.sip_pix2foc(wcs.wcs_world2pix(coo, 1), 1)
    weights = 1./events['response'] if response else None
    H, xedges, yedges = np.histogram2d(foc[:, 1]-0.5, foc[:, 0]-0.5, bins=imsz,
                                       range=([[0, imsz[0]], [0, imsz[1]]]),
                                       weights=weights)

    return H
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def integrate_map(band, skypos, tranges, skyrange, verbose=0, memlight=None,
                  hdu=None, retries=100, response=False, detsize=1.1):
    """
    Integrate an image over some number of time ranges. Use a reduced
	    memory optimization (at the expense of more web queries) if requested.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param tranges: Set of time ranges to retrieve the photon events, in
        GALEX time seconds.

    :type tranges: list

    :param skyrange: RA and Dec extents of the region of interest in degrees.

    :type skyrange: list

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param memlight: Reduce memory usage by breaking query into smaller
        segments of this size in seconds.

    :type memlight: float

    :param hdu: An existing HDU to modify.

    :type hdu: bool

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int

    :param response: Apply the response correction.

    :type response: bool

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :returns: numpy.ndarray - The integrated image.
    """

    imsz = gxt.deg2pix(skypos, skyrange)
    img = np.zeros(np.array(imsz, dtype='int32'))

    for trange in tranges:
        img += makemap(band, skypos, trange, skyrange,
            response=response, verbose=verbose, detsize=detsize)

    if response: # Intensity maps == average countrate maps.
        expt=np.sum([dbt.compute_exptime(band,trange) for trange in tranges])
        if expt>0:
            img/=expt

    return img
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def write_jpeg(filename, band, skypos, tranges, skyrange, stepsz=1.,
               overwrite=None, verbose=0, retries=100):
    """
    Write a 'preview' jpeg image from a count map.

    :param filename: Name of output jpeg to make.

    :type filename: str

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param tranges: Set of time ranges to retrieve the photon events,
        in GALEX time seconds.

    :type tranges: list

    :param skyrange: RA and Dec extent of the region of interest in degrees.

    :type skyrange: list

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
                                              verbose=verbose, retries=retries))

    return
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def movie(band, skypos, tranges, skyrange, framesz=0, verbose=0,
          memlight=None, coadd=False, response=False, hdu=None, retries=100,
          detsize=1.1):
    """
    Generate a movie (mov) file.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param tranges: Set of time ranges to retrieve the photon events,
        in GALEX time seconds.

    :type tranges: list

    :param skyrange: RA and Dec extents of the region of interest in degrees.

    :type skyrange: list

    :param framesz: The time bin size (depth) to use per frame, in seconds.

    :type framesz: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param memlight: Reduce memory usage by breaking query into smaller
        segments of this depth in seconds.

    :type memlight: float

    :param coadd: Integrated across all time ranges. (i.e. create a coadd)

    :type coadd: bool

    :param response: Apply the response correction.

    :type response: bool

    :param hdu: An existing HDU to modify.

    :type hdu: bool

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :returns: numpy.ndarray -- The movie file.
    """

    # Not defining stepsz creates a single full depth image.
    if verbose:
        print(tranges)

    if len(np.shape(tranges))==1:
        tranges=[tranges]

    if coadd or (len(tranges) == 1 and not framesz) or (not len(tranges)):
        if verbose > 2:
            print('Coadding across '+str(tranges))

        mv = integrate_map(band, skypos, tranges, skyrange,
                           verbose=verbose, memlight=memlight,
                           hdu=hdu, retries=retries, response=response,
                           detsize=detsize)
    else:
        mv = np.array(None) # Initialize to error gracefully in the case of no data
        for trange in tranges:
            stepsz = framesz if framesz else trange[1]-trange[0]
            try:
                steps = np.ceil((trange[1]-trange[0])/stepsz)
            except IndexError:
                return None # expt_raw == 0
            if framesz==0:
                trs=[trange]
            else:
                tsteps = np.arange(trange[0],trange[1],framesz)
                trs = np.vstack([tsteps,np.append(tsteps[1:],trange[1])]).T
            for tr in trs:
                t0,t1=tr
                img = integrate_map(band, skypos, [[t0, t1]], skyrange,
                                    verbose=verbose,
                                    memlight=memlight, hdu=hdu, retries=retries,
                                    response=response, detsize=detsize)
                if (img.min() == 0 and img.max() == 0) or (not
                                                np.isfinite(img).any()):
                    #if verbose > 1:
                    #    print('No data in frame {i}. Skipping...'.format(i=i))
                    continue
                else:
                    try:
                        mv.append(img)
                    except:
                        mv = [img]

    try:
        if len(np.where(mv>0)[0])==0:
            return np.array(None)
    except TypeError:
        pass

    return np.array(mv)

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def create_image(band, skypos, tranges, skyrange, framesz=0, verbose=0,
                 memlight=None, coadd=False, response=False, hdu=None,
                 retries=100, detsize=1.1):
    """
    Generate count or intensity images or movies at a given sky position and
        across given time ranges.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param tranges: Set of time ranges to retrieve the photon events,
        in GALEX time seconds.

    :type tranges: list

    :param skyrange: RA and Dec extents of the region of interest in degrees.

    :type skyrange: list

    :param framesz: The time bin size (depth) to use per frame, in seconds.

    :type framesz: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param memlight: Reduce memory usage by breaking query into smaller
        segments of this size in seconds.

    :type memlight: float

    :param coadd: Integrate across all time ranges. (i.e. make a coadd)

    :type coadd: bool

    :param response: Apply the relative response correction.

    :type response: bool

    :param hdu: An existing HDU to modify.

    :type hdu: bool

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :returns: numpy.ndarray -- The image file.
    """

    img = movie(band, skypos, tranges, skyrange, framesz=framesz,
                verbose=verbose, memlight=memlight,
                coadd=coadd, response=response, hdu=hdu, retries=retries,
                detsize=detsize)

    return np.array(img)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def write_images(band, skypos, tranges, skyrange, write_cnt=None,
                 write_int=None, framesz=0, verbose=0, memlight=None,
                 coadd=False, overwrite=None, retries=100,
                 write_cnt_coadd=False, write_int_coadd=False, detsize=1.1):
    """
    Generate a write various maps to files.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param tranges: Set of time ranges to retrieve the photon events,
        in GALEX time seconds.

    :type tranges: list

    :param skyrange: RA and Dec extents of the region of interest in degrees.

    :type skyrange: list

    :param write_cnt: Make count image?

    :type write_cnt: bool

    :param write_int: Make intensity image?

    :type write_int: bool

    :param framesz: The time bin size (depth) to use per frame, in seconds

    :type framesz: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param memlight: Reduce memory usage by breaking query into smaller
        segments of this size in seconds.

    :type memlight: float

    :param coadd: Integrate across all time ranges (i.e. create a coadd)

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

    for i in list(imtypes.keys()):
        if not imtypes[i]:
            continue

        img = create_image(band, skypos, tranges, skyrange, framesz=framesz,
            verbose=verbose, memlight=memlight, retries=retries, detsize=detsize,
            coadd=coadd or i in ['cnt_coadd','int_coadd'],
            response=i in ['int', 'int_coadd'])
        if img.tolist() is None:
            if verbose:
                print('No data found.')
            return

        # Add a conditional so that this is only created for multi-frame images
        tbl = movie_tbl(band, tranges, framesz=framesz, verbose=verbose,
                        retries=retries) if i in ['int', 'int_coadd'] else False

        hdu = pyfits.PrimaryHDU(img)
        hdu = fits_header(band, skypos, tranges, skyrange, verbose=verbose,
                          hdu=hdu, retries=retries)

        hdulist = pyfits.HDUList([hdu, tbl]) if tbl else pyfits.HDUList([hdu])

        if verbose:
            print('Writing image to {o}'.format(o=imtypes[i]))

        hdulist.writeto(imtypes[i], overwrite=overwrite)

    return
# ------------------------------------------------------------------------------
