"""
.. module:: MCUtils
   :synopsis: A library of generic utilities programs that C. Million wants to
       keep separate in order to re-use across other projects.

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
from astropy.io import fits as pyfits
from builtins import str
import numpy as np
import requests
from sys import stdout

# ------------------------------------------------------------------------------
def area(radius):
    """
    Returns the area of a circle with a given radius.

    :param radius: The radius of the cricle.

    :type radius: float

    :returns: float -- The area of the circle.
    """

    return np.pi*radius**2.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def distance(a, b, c, d):
    """
    Computes Euclidean distance between [a,b] and [c,d].

    :param a: x-coordinate of first data point.

    :type a: float

    :param b: y-coordinate of first data point.

    :type b: float

    :param c: x-coordinate of second data point.

    :type c: float

    :param d: y-coordinate of second data point.

    :type d: float

    :returns: float -- The Euclidean distance between the two points.
    """

    return np.sqrt((a-c)**2. + (b-d)**2.)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def rotvec(vector, theta):
    """
    Rotate vectors clockwise by theta degrees.

    :param vector: The vector to rotate. Must have (2,n) shape.

    :type vector: numpy.ndarray

    :param theta: Angle to rotate the vector, in degrees. Must have (n,) shape.

    :type theta: numpy.ndarray

    :returns: numpy.ndarray -- The rotated vector with (2,n) shape.
    """

    R = np.array([[np.cos(np.radians(theta)), -np.sin(np.radians(theta))],
                  [np.sin(np.radians(theta)), np.cos(np.radians(theta))]])

    return np.dot(R.T, vector)[0]
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def rms(data):
    """
    Return the root-mean-square of the set of values.

    :param data: The set of values from which to calculate the
        root-mean-squared.

    :type data: list

    :returns: float -- The root-mean-square of the set of values.
    """

    return np.sqrt(np.square(np.array(data)-np.array(data).mean()).mean())
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def print_inline(text, blanks=60):
    """
    For updating text in place without a carriage return.

    :param text: Message to print to standard out.

    :type text: str

    :param blanks: Number of white spaces to prepend to message.

    :type blanks: int
    """

    stdout.write(" "*blanks+"\r")
    stdout.write(str(str(text)+'\r'))
    stdout.flush()
    return
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def manage_requests2(query, maxcnt=100, wait=1, timeout=60., verbose=0):
    """
    Make simple 'requests' calls more robust against network issues.

    :param query: The URL containing the query.

    :type query: str

    :param maxcnt: The maximum number of attempts to make before failure.

    :type maxcnt: int

    :param wait: The length of time to wait before attempting the query again.
        Currently a placeholder.

    :type wait: int

    :param timeout: The length of time to wait for the server to send data
        before giving up, specified in seconds.

    :type timeout: float

    :param verbose: If > 0, print additional messages to STDOUT. Higher value
        represents more verbosity.

    :type verbose: int

    :returns: requests.Response or None -- The response from the server. If the
        query does not receive a response, returns None.
    """

    # Keep track of the number of failures.
    cnt = 0

    # This will keep track of whether we've gotten at least one
    # successful response.
    successful_response = False

    while cnt < maxcnt:
        try:
            r = requests.get(query, timeout=timeout)
            successful_response = True
        except requests.exceptions.ConnectTimeout:
            if verbose:
                print("Connection time out.")
            cnt += 1
            continue
        except requests.exceptions.ConnectionError:
            if verbose:
                print("Domain does not resolve.")
            cnt += 1
            continue
        except:
            if verbose:
                print('bad query? {q}'.format(q=query))
            cnt += 1
            continue
        if r.json()['status'] == 'EXECUTING':
            if verbose > 1:
                print_inline('EXECUTING')
            cnt = 0
            continue
        elif r.json()['status'] == 'COMPLETE':
            if verbose > 1:
                print_inline('COMPLETE')
            break
        elif r.json()['status'] == 'ERROR':
            print('ERROR')
            print('Unsuccessful query: {q}'.format(q=query))
            raise ValueError(r.json()['msg'])
        else:
            print('Unknown return: {s}'.format(s=r.json()['status']))
            cnt += 1
            continue

    if not successful_response:
        # Initiate an empty response object in case
        # the try statement is never executed.
        # r = requests.Response()
        r = None

    return r
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def find_nearest_lower(array, value):
    """
    Finds the index of the value in the array that is closest without going
        over. This method assumes that:
            1. 'value' is within the range of 'array'.
            2. 'array' is ordered.
            3. 'array' has no gaps.

    :param array: Array of values to search.

    :type array: numpy.ndarray

    :param value: Value to find the closest match without going over.

    :type value: float

    :returns: int -- The index of the array element closest to value without
        going over.
    """

    idx = (np.abs(array-value)).argmin()

    if array[idx] > value:
        idx -= 1

    return idx
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def get_fits_data(filename, dim=0, verbose=0):
    """
    Reads FITS data. A wrapper for common pyfits commands.

    :param filename: The name of the FITS file to retrieve the data from.

    :type filename: str

    :param dim: The extension to retrieve the data from, 0=Primary, 1=First
        Extension, etc.

    :type dim: int

    :param verbose: If > 0, print messages to STDOUT.

    :type verbose: int

    :returns: Data instance -- The data from the 'dim' HDU.
    """

    if verbose:
        print("         ", filename)

    hdulist = pyfits.open(filename, memmap=1)

    data = hdulist[dim].data

    hdulist.close()

    return data
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def get_fits_header(filename):
    """
    Reads a FITS header. A wrapper for common astropy.io.fits commands.

    :param filename: The name of the FITS file to retrieve the header from.

    :type filename: str

    :returns: Header instance -- The header from the primary HDU.
    """

    hdulist = pyfits.open(filename, memmap=1)

    htab = hdulist[0].header

    hdulist.close()

    return htab
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def get_tbl_data(filename, comment='|'):
    """
    Reads data from a table into a numpy array.

    :param filename: The name of the FITS file to read.

    :type filename: str

    :param comment: The symbol that represents a comment.

    :type comment: str

    :returns: numpy.ndarray -- The table data.
    """

    f = open(filename)
    lines = f.readlines()
    tbl = []

    for line in lines:
        if line[0] != comment:
            strarr = str.split(str(line))
            if len(strarr) > 0:
                tbl.append(strarr)

    return np.array(tbl, dtype='float64')
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def angularSeparation(ra1, dec1, ra2, dec2):
    """
    Compute angular separation in degrees of points on the sky.
        It is important, especially for small angular separations, that the
        values for ra[01] and dec[01] have precision of float64 or better.
        Now uses the haversine formula which is stable for small angles.

    :param ra1: The right ascension of the first coordinate.

    :type ra1: float

    :param dec1: The declination of the first coordinate.

    :type dec1: float

    :param ra2: The right ascension of the second coordinate.

    :type ra2: float

    :param dec2: The declination of the second coordinate.

    :type dec2: float

    :returns: float -- The angular separation, in degrees, on the sky.
    """

    d2r = np.pi/180.
    ra2deg = 1./d2r

    d1 = dec1*d2r
    d2 = dec2*d2r

    r1 = ra1*d2r
    r2 = ra2*d2r

    a = np.sin((d2-d1)/2.)**2.+np.cos(d1)*np.cos(d2)*np.sin((r2-r1)/2.)**2.
    r = 2*np.arcsin(np.sqrt(a))

    return r*ra2deg
# ------------------------------------------------------------------------------
