"""
.. module:: PhotonTools

   :synopsis: @CHASE - Please write the synopsis for this module.@

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

# @CHASE - We need to avoid "import *", can we do specific imports here?@
import csv
from astropy.io import fits as pyfits
import numpy as np
import math
from astropy import wcs as pywcs
from FileUtils import *
from CalibrationTools import *
from CalUtils import *
from MCUtils import *

# ------------------------------------------------------------------------------
def scanraw6(raw6file):
    """
    Test of raw6 files.
    @CHASE - Notes suggest this is not used anywhere, can remove?@

    :param raw6file: The name of the raw6 file to read.

    :type raw6file: str
    """

    raw6hdulist = pyfits.open(raw6file, memmap=1)
    raw6htab = raw6hdulist[1].header
    nphots = raw6htab['NAXIS2']
    print "			nphots= ", nphots

    for i in xrange(nphots):
        t = np.array(raw6hdulist[1].data.field('t')[i])
        phb1 = np.array(raw6hdulist[1].data.field('phb1')[i], dtype='int64')
        phb2 = np.array(raw6hdulist[1].data.field('phb2')[i], dtype='int64')
        phb3 = np.array(raw6hdulist[1].data.field('phb3')[i], dtype='int64')
        phb4 = np.array(raw6hdulist[1].data.field('phb4')[i], dtype='int64')
        phb5 = np.array(raw6hdulist[1].data.field('phb5')[i], dtype='int64')

        if t == 874219775.53 or t == 874219893.195 or t == 874220274.02:
            print t, phb1, phb2, phb3, phb4, phb5

    return
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def read_xfile(xfile):
    """
    Reads an 'x-file', which is what the GALEX mission team called their
    photon list files.

    :param xfile: Name of the x-file to read.

    :type xfile: str

    :returns: tuple -- A three-element tuple containing the timestamps, RA, and
    DEC values in the 'x-file'.
    """

    print "Reading xfile ", xfile
    hdulist = pyfits.open(xfile, memmap=1)

    scale = 2147483647.
    flag = np.array(hdulist[1].data.field('flags'))
    t = np.array(hdulist[1].data.field('t'))
    X = np.array(hdulist[1].data.field('X'))/scale
    Y = np.array(hdulist[1].data.field('Y'))/scale
    Z = np.array(hdulist[1].data.field('Z'))/scale
    skip = -2147483647
    ix = ((X > skip)).nonzero()

    dec = np.arcsin(Z)
    ra = np.arctan2(Y, X)

    return t, ra, dec
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def get_ra_dec(fitsfile):
    """
    Pulls the RA_CENT and DEC_CENT from a FITS header.

    :param fitsfile: The name of the FITS file to read.

    :type fitsfile: str

    :returns: tuple -- A two-element tuple containing the RA and DEC.
    """

    print "Reading boresite center RA, Dec from ", fitsfile
    hdulist = pyfits.open(fitsfile, memmap=1)
    ra = hdulist[0].header['RA_CENT']
    dec = hdulist[0].header['DEC_CENT']

    print '		', ra, ', ', dec

    hdulist.close()

    return ra, dec
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def readfromcsv(csvfile):
    """
    Reads out t, ra, dec from a photon CSV file.
    @CHASE - Question of whether this is redundant/unnecessary, shall we
    remove?@

    :param csvfile: Name of the CSV file to read.

    :type csvfile: str

    :returns: list -- A list containing the timestamps, RA, and DEC from the
    CSV.
    """

    reader = csv.reader(open(csvfile, 'rb'), delimiter=',', quotechar='|')

    t, ra, dec = [], [], []

    for row in reader:
        t.append(float(row[0]))
        ra.append(float(row[8]))
        dec.append(float(row[9]))

    return [t, ra, dec]
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def makeimage(t, ra, dec, stepsz, dims, racent, deccent):
    """
    Generates a count (cnt) image based on t, ra, dec values.

    :param t: Set of timestamps.

    :type t: list @CHASE - please confirm data type.@

    :param ra: Set of right ascension values.

    :type ra: list @CHASE - please confirm data type.@

    :param dec: Set of declination values.

    :type dec: list @CHASE - please confirm data type.@

    :param stepsz: The desired step size, in seconds, for each image frame.
    @CHASE - Please check this description, is it per frame, or some sort of
    chunk size as the image is built, etc.?@

    :type stepsz: int @CHASE - can this be a float?@

    :param dims: @CHASE - please provide description for this parameter.@

    :type dims: @CHASE - please specify data type.@

    :param racent: The central RA of the image, in degrees. @CHASE - please
    confirm.@

    :type racent: float

    :param deccent: The central DEC of the image, in degrees. @CHASE - please
    confirm.@

    :type deccent: float

    :returns: list -- A list of 2D images.
    """

    t, ra, dec = np.array(t), np.array(ra), np.array(dec)
    mint, maxt = min(t), max(t)
    trange = maxt-mint
    if stepsz > (maxt-mint):
        stepsz = 0

    if stepsz <= 0 or stepsz >= trange:
        image, xedges, yedges = np.histogram2d(
            dec, ra, bins=dims, range=([deccent-0.8, deccent+0.8],
                                       [racent-0.8, racent+0.8]))
        return np.fliplr(image)
    else:
        if not math.fmod(trange, stepsz):
            steps = math.floor(trange/stepsz)+1
        else:
            steps = trange/stepsz

    ims = []
    for i in xrange(steps):
        if (mint+(i+1)*stepsz) < maxt:
            ix = ((t >= (mint+i*stepsz)) & (t < (mint+(i+1)*stepsz))).nonzero()
        else:
            ix = (t >= (mint+i*stepsz)).nonzero()

        # Fix the bin edges on this so that it can't wobble
        image, xedges, yedges = np.histogram2d(dec[ix], ra[ix], bins=dims,
                                               range=([deccent-0.8,
                                                       deccent+0.8],
                                                      [racent-0.8,
                                                       racent+0.8]))
        image = np.fliplr(image)
        # Numpy appears to make it stupidly difficult to stack images in a way
        # that would be straightforward in IDL... hence the ridiculous logic
        # structure below.
        if i == 0:
            ims = image
        elif i == 1:
            ims = np.concatenate(([ims], [image]))
        else:
            ims = np.concatenate((ims, [image]))

    return ims
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def readrangesfromcsv(csvfile):
    """
    Read a CSV file and return the minimum and maximum time, right ascension,
    and declination. Also returns the total number of rows.

    :param csvfile: The CSV file to read.

    :type csvfile: str

    :returns: tuple -- A 7-element tuple containing the min. timestamp, max.
    timestamp, min. RA, max. RA, min. DEC, max. DEC, and total number of rows.
    """

    reader = csv.reader(open(csvfile, 'rb'), delimiter=',', quotechar='|')

    t, ra, dec = [], [], []

    for row in reader:
        if row[0] and row[8] and row[9] and not float(row[10]):
            t.append(float(row[0]))
            ra.append(float(row[8]))
            dec.append(float(row[9]))

    return min(t), max(t), min(ra), max(ra), min(dec), max(dec), len(t)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def makeimage_memlight(csvfile, stepsz, dims, racent, deccent):
    """
    A memory light way of making count (cnt) from a photon CSV file. This is a
    slower but more memory efficient use of 'makeimage'. Try this if you run out
    of memory while using 'makeimage' (which could well happen on NUV datasets).
    The stepsz (e.g. movie time resolution) is not yet implemented.

    :param csvfile: Name of the CSV file to use.

    :type csvfile: str

    :param stepsz: The desired step size, in seconds, for each image frame.
    @CHASE - Please check this description, also, since it's not implemented,
    should this just be removed for now and put in when/if it's implemented?@

    :type stepsz: int @CHASE - can this be a float?@

    :param dims: @CHASE - please provide description for this parameter, looks
    like the size of the image (in pixels?).@

    :type dims: @CHASE - please specify data type.@

    :param racent: The central RA of the image, in degrees. @CHASE - please
    confirm.@

    :type racent: float

    :param deccent: The central DEC of the image, in degrees. @CHASE - please
    confirm.@

    :type deccent: float

    :returns: numpy.ndarray -- The image as a 2D array.
    """

    print "Getting space and time ranges..."
    mint, maxt, minra, maxra, mindec, maxdec, rowcount = readrangesfromcsv(
        csvfile)
    print "		time: ", mint, " : ", maxt
    print "		  ra: ", minra, " : ", maxra
    print "		 dec: ", mindec, " : ", maxdec

    print "Allocating memory to image array..."
    image = np.zeros([dims, dims], dtype=float32)

    print "Binning photon data..."
    reader = csv.reader(open(csvfile, 'rb'), delimiter=',', quotechar='|')
    chunksz = 100000
    cnt, chunk, n = 0, 0, 0
    t, ra, dec = [], [], []
    for row in reader:
        n += 1
        if cnt < chunksz:
            cnt += 1
            t.append(float(row[0]))
            ra.append(float(row[8]))
            dec.append(float(row[9]))

        if cnt >= chunksz or n == rowcount:
            chunk += 1
            print_inline("Coadding chunk "+str(chunk))
            hist, xedges, yedges = np.histogram2d(
                dec, ra, bins=dims, range=([deccent-0.8, deccent+0.8],
                                           [racent-0.8, racent+0.8]))
            image += hist
            cnt = 0
            t, ra, dec = [], [], []

    print n

    return np.fliplr(image)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def writephots2img(image, header, racent, deccent):
    """
    @CHASE - This method should just be removed entirely.@
    """

    pyfits.writeto('exp.fits', image, header)

    return 0
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def createheader(rarange, decrange, npix, racent, deccent, expstart=0.,
                 expend=0., exptime=0., header=pyfits.Header()):
    """
    Populates a FITS header for a GALEX image.

    :param rarange: @CHASE - please describe.@

    :type rarange: @CHASE - float? list?@

    :param decrange: @CHASE - please describe.@

    :type decrange: @CHASE - float? list?@

    :param npix: Number of pixels in one dimension @CHASE - does this assume
    a square image?@

    :type npix: int

    :param racent: RA center of the image, in degrees.

    :type racent: float

    :param deccent: DEC center of the image, in degrees.

    :type deccent: float

    :param expstart: Start time of the exposure @CHASE - in JD? GALEX time?@

    :type expstart: float @CHASE - int/long/float?@

    :param expend: End time of the exposure @CHASE - in JD? GALEX time?@

    :type expend: float @CHASE - int/long/float?@

    :param exptime: Total exposure time, in seconds.

    :type exptime: float

    :param header: Header object to populate.

    :type header: astropy.io.fits.Header object.

    :returns: astropy.io.fits.Header -- The updated Header object.
    """

    # The following values will be set automatically by pyfits
    # header.update(key='simple', value='T',
    #   comment='conforms to FITS standard')
    # header.update(key='bitpix', value=-64, comment='array data type')
    # The following values will be set automatically by pyfits
    # header.update(key='naxis', value=2, comment='number of array dimensions')
    # header.update(key='naxis1', value=1000)
    # header.update(key='naxis2', value=1000)
    # header.update(key='extend', value='T')

    header.update(key='extend', value='T',
                  comment='FITS dataset may contain extensions')
    header.update(key='cdelt1', value=rarange/npix)
    header.update(key='cdelt2', value=decrange/npix)
    header.update(key='equinox', value=2000.)
    header.update(key='epoch', value=2000.)
    header.update(key='ctype1', value='RA---TAN')
    header.update(key='ctype2', value='DEC--TAN')
    header.update(key='crpix1', value=(npix/2.)+.5)
    header.update(key='crpix2', value=(npix/2.)+.5)
    header.update(key='crval1', value=racent)
    header.update(key='crval2', value=deccent)
    header.update(key='crota2', value=0.)
    header.update(key='bunit', value='        ')
    header.update(key='bscale', value=1.)
    header.update(key='bzero', value=0.)
    header.update(key='ra_cent', value=racent)
    header.update(key='dec_cent', value=deccent)
    header.update(key='twist', value=-999.)
    header.update(key='expstart', value=expstart)
    header.update(key='expend', value=expend)
    if not exptime:
        header.update(key='exptime', value=expend-expstart)
    else:
        header.update(key='exptime', value=exptime)

    return header
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def create_cnt(csvfile, imsz, cntfile, outfile, matchtimes=0, expstart=0,
               expend=0):
    """
    Writes a count (cnt) image to a FITS file from a photon CSV file.

    :param csvfile: Name of photon event CSV file to use.

    :type csvfile: str

    :param imsz: Size of the image @CHASE - looks like it must be rectangular?@

    :type imsz: int

    :param cntfile: Name of count image file to use.

    :type cntfile: str

    :param outfile: Name of output file to create.

    :type outfile: str

    :param matchtimes: @CHASE - please describe this parameter.@

    :type matchtimes: int @CHASE - This should be a Boolean (True/False)@

    :param expstart: Start time of the image @CHASE - JD? GALEX time?@

    :type expstart: @CHASE - int/long/float?@

    :param expend: End time of the image @CHASE - JD? GALEX time?@

    :type expend: @CHASE - int/long/float?@
    """

    hdulist = pyfits.open(cntfile)
    hdr = hdulist[0].header
    hdulist.close()

    wcs = pywcs.WCS(naxis=2)
    wcs.wcs.crpix = [hdr['crpix1'], hdr['crpix2']]
    wcs.wcs.cdelt = [hdr['cdelt1'], hdr['cdelt2']]
    wcs.wcs.crval = [hdr['crval1'], hdr['crval2']]
    wcs.wcs.ctype = [hdr['ctype1'], hdr['ctype2']]
    racent, deccent = hdr['crval1'], hdr['crval2']
    if not expend and not expstart:
        expstart, expend = hdr['expstart'], hdr['expend']
    elif expend and not expstart:
        expstart = hdr['expstart']
    elif expstart and not expend:
        expend = hdr['expend']
    exptime = hdr['exptime']

    print 'Time range of ['+str(expstart)+', '+str(expend)+'].'

    if hdr['band'] == 1:
        band = 'NUV'
    elif hdr['band'] == 2:
        band = 'FUV'
    else:
        print "Band not specified in header."
    eclipse = hdr['eclipse']

    print wcs.wcs.name

    image = np.zeros([imsz, imsz], dtype=np.dtype('f'))

    reader = csv.reader(open(csvfile, 'rb'), delimiter=',', quotechar='|')

    mint, maxt = 0, 0
    chunksz = 100000
    cnt, chunk, n = 0, 0, 0
    t, coo = [], []
    for row in reader:
        n += 1
        if cnt < chunksz:
            cnt += 1
            if row[0] and row[8] and row[9] and not float(row[10]):
                if (matchtimes and ((float(row[0])/1000.)+GPSSECS < expstart or
                                    (float(row[0])/1000.)+GPSSECS > expend)):
                    continue
                coo.append([np.float64(row[8]), np.float64(row[9])])
                t.append(np.float64(row[0])/1000)

        if cnt >= chunksz:
            chunk += 1
            if len(t) > 0:
                pix = wcs.wcs_world2pix(coo, 1)
                foc = wcs.sip_pix2foc(pix, 1)
                xpix = foc[:, 0]
                ypix = foc[:, 1]
                print_inline("Coadding chunk "+str(chunk))
                H, xedges, yedges = np.histogram2d(ypix-0.5, xpix-0.5,
                                                   bins=imsz, range=([0, imsz],
                                                                     [0, imsz]))
                image += H

                # Find exposure time range
                if min(t) < mint or mint == 0:
                    mint = min(t)
                if max(t) > maxt:
                    maxt = max(t)

            else:
                print_inline("Chunk "+str(chunk)+" contained no unflagged data.")
            cnt = 0
            t, coo = [], []

    if t:
        pix = wcs.wcs_world2pix(coo, 1)
        foc = wcs.sip_pix2foc(pix, 1)
        xpix = foc[:, 0]
        ypix = foc[:, 1]
        if xpix.any() and ypix.any():
            print_inline("Coadding final chunk "+str(chunk+1))
            H, xedges, yedges = np.histogram2d(ypix-0.5, xpix-0.5, bins=imsz,
                                               range=([0, imsz], [0, imsz]))
            image += H
        else:
            print "Chunk ", chunk, " contained no unflagged data."

    print "Finished reading ", n, " photons."

    header = createheader(-1.6, 1.6, imsz, racent, deccent,
                          expstart=mint+GPSSECS, expend=maxt+GPSSECS,
                          exptime=exptime, header=hdr)
    pyfits.writeto(outfile, image, header, clobber=True)

    return
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def create_dose(csvfile, imsz, cntfile, outfile):
    """
    Writes a dose map (detector space image) from a photon CSV file.

    :param csvfile: Name of photon event CSV file to use.

    :type csvfile: str

    :param imsz: Size of the image @CHASE - looks like it must be rectangular?@

    :type imsz: int

    :param cntfile: Name of count image file to use.

    :type cntfile: str

    :param outfile: Name of output file to create.

    :type outfile: str
    """

    hdulist = pyfits.open(cntfile)
    hdr = hdulist[0].header
    hdulist.close()

    print "Getting space and time ranges..."
    mint, maxt, minra, maxra, mindec, maxdec, rowcount = readrangesfromcsv(
        csvfile)
    print "		time: ", mint, " : ", maxt
    print "		  ra: ", minra, " : ", maxra
    print "		 dec: ", mindec, " : ", maxdec

    image = np.zeros([imsz, imsz], dtype=np.dtype('i'))

    reader = csv.reader(open(csvfile, 'rb'), delimiter=',', quotechar='|')

    chunksz = 100000
    cnt, chunk, n = 0, 0, 0
    t, xi, eta = [], [], []
    for row in reader:
        n += 1
        if cnt < chunksz:
            cnt += 1
            t.append(np.float64(row[0]))
            xi.append(np.float64(row[6]))
            eta.append(np.float64(row[7]))

        if cnt >= chunksz:
            chunk += 1
            print_inline("Coadding chunk "+str(chunk))
            H, xedges, yedges = np.histogram2d(xi, eta, bins=imsz,
                                               range=([-25000, 25000],
                                                      [-25000, 25000]))
            image += H
            t, xi, eta = xi, eta

    if t:
        print_inline("Coadding final chunk "+str(chunk+1))
        H, xedges, yedges = np.histogram2d(xi, eta, bins=imsz,
                                           range=([-25000, 25000],
                                                  [-25000, 25000]))
        image += H

    print "Finished reading ", n, " photons."

    racent, deccent = hdr['crval1'], hdr['crval2']
    header = createheader(-1.6, 1.6, imsz, racent, deccent)
    pyfits.writeto(outfile, H, header, clobber=True)

    # @CHASE - Can this just return nothing?@
    return 0
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def load_stims(csvfile, band, eclipse):
    """
    Returns only the stim events from a photon CSV file.

    :param csvfile: Name of photon event CSV file to use.

    :type csvfile: str

    :param band: The band to load stims for, either 'FUV' or 'NUV'.

    :type band: str

    :param eclipse: The eclipse number to load stims for.

    :type eclipse: int @CHASE - Confirm this is int/long/float.@

    :returns: tuple -- A four-element tuple containing arrays of the times,
    detector x, detector y, and stim ID (1, 2, 3 or 4) of the likely stim
    events.
    """

    reader = csv.reader(open(csvfile, 'rb'), delimiter=',', quotechar='|')

    x, y, t = [], [], []

    for row in reader:
        x.append(np.float64(row[1]))
        y.append(np.float64(row[2]))

    return find_stims(x, y, band, eclipse)
# ------------------------------------------------------------------------------
