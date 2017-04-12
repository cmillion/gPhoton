"""
.. module:: FileUtils
   :synopsis: Methods for reading data from files or web queries.

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
from astropy.io import fits as pyfits
from builtins import str
from builtins import range
import os
import numpy as np
# gPhoton imports.
import gPhoton.cal as cal
import gPhoton.gQuery as gQuery
from gPhoton.MCUtils import manage_requests2

def get_raw_paths(eclipse,verbose=0):
    url = gQuery.raw_data_paths(eclipse)
    if verbose>1:
        print(url)
    r = manage_requests2(url)
    out = {'NUV':None,'FUV':None,'scst':None}
    for f in r.json()['data']['Tables'][0]['Rows']:
        if (f[1].strip()=='NUV') or (f[1].strip()=='FUV'):
            out[f[1].strip()]=f[2]
        elif f[1].strip()=='BOTH': # misnamed scst path
            out['scst']=f[2]
    return out

def download_data(eclipse,band,ftype,datadir='./',verbose=0):
    urls = get_raw_paths(eclipse,verbose=verbose)
    if not ftype in ['raw6','scst']:
        raise ValueError('ftype must be either raw6 or scst')
    if not band in ['NUV','FUV']:
        raise ValueError('band must be either NUV or FUV')
    url = urls[band] if ftype is 'raw6' else urls['scst']
    if not url:
        print('Unable to locate {f} file on MAST server.'.format(f=ftype))
        return None
    if url[-6:]=='.gz.gz': # Handling a very rare mislabeling of the URL.
        url = url[:-3]
    if verbose>1:
        print(url)
    if not datadir:
        datadir = '.'
    if datadir and datadir[-1]!='/':
        datadir+='/'
    filename = url.split('/')[-1]
    opath = '{d}{f}'.format(d=datadir,f=filename)
    if os.path.isfile(opath):
        print('Using {f} file already at {d}'.format(
                                    f=ftype,d=os.path.abspath(opath)))
    else:
        # FIXME: All exceptions are treated the same. Could be missing data or
        #  a network outage. Should treat these cases differently.
        try:
            cal.download_with_progress_bar(url,opath)
        except:
            print('Unable to download data from {u}'.format(u=url))
            raise
            #return None
    return opath


# ------------------------------------------------------------------------------
def load_raw6(raw6file):
    """
    Reads a raw6 file. Just wraps some pyfits I/O commands.

    :param raw6file: Name of the raw6 file to read.

    :type raw6file: str

    :returns: tuple -- A two-element tuple containing the header of the
        first extension, and the HDUList object returned from astropy.io.fits.
    """

    print("		", raw6file)
    hdulist = pyfits.open(raw6file, memmap=1)
    htab = hdulist[1].header
    hdulist.close()
    return htab, hdulist
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def load_aspect(aspfile):
    """
    Loads a set of aspect files into a bunch of arrays.

    :param aspfile: List of aspect files (+paths) to read.

    :type aspfile: list

    :returns: tuple -- Returns a six-element tuple containing the RA, DEC,
        twist (roll), time, header, and aspect flags. Each of these EXCEPT for
        header, are returned as numpy.ndarrays. The header is returned as a dict
        containing the RA, DEC, and roll from the headers of the aspec files in
        numpy.ndarrays.
    """

    ra, dec = np.array([]), np.array([])
    twist, time, aspflags = np.array([]), np.array([]), np.array([])

    header = {'RA':[], 'DEC':[], 'ROLL':[]}
    for i in range(len(aspfile)):
        print("         ", aspfile[i])
        hdulist = pyfits.open(aspfile[i], memmap=1)
        ra = np.append(ra, np.array(hdulist[1].data.field('ra')))
        dec = np.append(dec, np.array(hdulist[1].data.field('dec')))
        twist = np.append(twist, np.array(hdulist[1].data.field('roll')))
        time = np.append(time, np.array(hdulist[1].data.field('t')))
        aspflags = np.append(aspflags,
                             np.array(hdulist[1].data.field('status_flag')))
        header['RA'] = np.append(header['RA'],
                                 np.zeros(len(hdulist[1].data.field('ra')))+
                                 hdulist[0].header['RA_CENT'])
        header['DEC'] = np.append(header['DEC'],
                                  np.zeros(len(hdulist[1].data.field('dec')))+
                                  hdulist[0].header['DEC_CENT'])
        header['ROLL'] = np.append(header['ROLL'],
                                   np.zeros(len(hdulist[1].data.field('roll')))+
                                   hdulist[0].header['ROLL'])

    # Sort by time.
    ix = np.argsort(time)

    hdulist.close()

    return ra[ix], dec[ix], twist[ix], time[ix], header, aspflags[ix]
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def web_query_aspect(eclipse, retries=20):
    """
    Grabs the aspect data from MAST databases based on eclipse.

    :param eclipse: The number of the eclipse to retrieve aspect files for.

    :type eclipse: int

    :param retries: The number of times to retry a query before giving up.

    :type retries: int

    :returns: tuple -- Returns a six-element tuple containing the RA, DEC,
        twist (roll), time, header, and aspect flags. Each of these EXCEPT for
        header, are returned as numpy.ndarrays. The header is returned as a dict
        containing the RA, DEC, and roll from the headers of the aspec files in
        numpy.ndarrays.
    """

    print("Attempting to query MAST database for aspect records.")
    entries = gQuery.getArray(gQuery.aspect_ecl(eclipse), retries=retries)
    n = len(entries)
    print('		Located '+str(n)+' aspect entries.')
    if not n:
        print("No aspect entries for eclipse "+str(eclipse))
        return
    ra, dec, twist, time, flags = [], [], [], [], []
    header = {'RA':[], 'DEC':[], 'ROLL':[]}
    ra0, dec0, twist0 = [], [], []
    for i in range(n):
        # The times are *1000 in the database to integerify
        time.append(float(entries[i][2])/1000.)
        ra.append(float(entries[i][3]))
        dec.append(float(entries[i][4]))
        twist.append(float(entries[i][5]))
        flags.append(float(entries[i][6]))
        ra0.append(float(entries[i][7]))
        dec0.append(float(entries[i][8]))
        twist0.append(float(entries[i][9]))

    # Need to sort the output so that it is time ordered before returning.
    # Although it should already be ordered by time because that is requested
    #  in the SQL query above. If this is time consuming, remove it.
    ix = np.argsort(np.array(time))
    header = ({'RA':np.array(ra0)[ix], 'DEC':np.array(dec0)[ix],
               'ROLL':np.array(twist0)[ix]})

    return (np.array(ra)[ix], np.array(dec)[ix], np.array(twist)[ix],
            np.array(time)[ix], header, np.array(flags)[ix])
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def create_ssd_filename(band, eclipse):
    """
    Returns the Stim Separation Data (SSD) calibration file name.

    :param band: The band to create the SSD for, either 'FUV' or 'NUV'.

    :type band: str

    :param eclipse: The eclipse number to create the SSD file for.

    :type eclipse: int

    :returns: str -- The name of the SSD file to create.
    """

    return "SSD_"+band.lower()+"_"+str(eclipse)+".tbl"
# ------------------------------------------------------------------------------
