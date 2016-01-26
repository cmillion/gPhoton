"""
.. module:: MCUtils

   :synopsis: A library of generic utilities programs that C. Million wants to
   keep separate in order to re-use across other projects.

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

from sys import stdout
import numpy as np
import requests
import time
from astropy.io import fits as pyfits
import csv

# ------------------------------------------------------------------------------
def error(program, note=''):
    """
    Print ERROR with callback and note.

    :param program: The method associated with the error.
    @CHASE - please confirm.@

    :type program: function @CHASE - not sure how to describe this data type.@

    :param note: Additional text to include with error message.

    :type note: str
    """

    print '*** ERROR: '+str(program)

    if note:
        print '    '+str(note)

    return
# ------------------------------------------------------------------------------

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

    :param data: The set of values to calculate the root-mean-square.

    :type data: list @CHASE - confirm data type.@

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
def manage_requests2(query, maxcnt=100, wait=10, timeout=60., verbose=0):
    """
    Make simple 'requests' calls more robust against network issues.

    :param query: The URL containing the query.

    :type query: str

    :param maxcnt: The maximum number of attempts to make before failure.

    :type maxcnt: int

    :param wait: The length of time to wait before attempting the query again.
    @CHASE - This parameter is not used in manage_requests2, remove?@

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
        except requests.exceptions.ConnectionError as e:
            if verbose:
                print "Domain does not resolve."
            cnt += 1
            continue
        except requests.exceptions.ConnectTimeout as e:
            if verbose:
                print "Connection time out."
            cnt += 1
            continue
        except:
            if verbose:
                print 'bad query? {q}'.format(q=query)
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
            print 'ERROR'
            print 'Unsuccessful query: {q}'.format(q=query)
            raise ValueError(r.json()['msg'])
        else:
            print 'Unknown return: {s}'.format(s=r.json()['status'])
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
def manage_requests(query, maxcnt=100, wait=10, timeout=60.):
    """
    Make simple 'requests' calls more robust against network issues.
    @CHASE - Is it still worth keeping this old method around?@

    :param query: The URL containing the query.

    :type query: str

    :param maxcnt: The maximum number of attempts to make before failure.

    :type maxcnt: int

    :param wait: The length of time to wait before attempting the query again.

    :type wait: int

    :param timeout: The length of time to wait for the server to send data
    before giving up, specified in seconds.

    :type timeout: float
    """

    # Keep track of the number of failures.
    cnt = 0

    while cnt < maxcnt:
        try:
            r = requests.get(query, timeout=timeout)
            # NOTE: This specifically tests for the return structure
            # from the MAST photon database and is therefore
            # not a good general test condition.
            test = r.json()['Tables'][0]['Rows']
            return r
        except:
            # This except, which doesn't raise, gets a pass because I really
            # do want it to catch every possible exception, and it will raise
            # eventually below with the unsuccessful query print statements.
            time.sleep(wait)
            cnt += 1
            print_inline('Query retry attempt '+str(int(cnt))+'.')

    print 'Query unsuccessful after '+str(int(maxcnt))+' attempts.'
    print '		'+str(query)

    return False
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def wheretrue(conditions):
    """
    Returns indices for which the input conditions are true. Example:

    .. code-block:: python

        a = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
	    wheretrue(a > 5)
	    array([5, 6, 7, 8])
	    a[wheretrue(a > 5)]
	    array([6, 7, 8, 9])

    :param conditions: Expression that must evaluate to True/False.

    :type conditions: numpy.ndarray

    :returns: numpy.ndarray -- Array of indexes where the condition is True.
    """

    return np.where(conditions == True)[0]
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def wherefalse(conditions):
    """
    Returns indices for which the input conditions are false.

    .. code-block:: python

	    a = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
	    wherefalse(a > 5)
	    array([0, 1, 2, 3, 4])
	    a[wherefalse(a > 5)]
	    array([1, 2, 3, 4, 5])

    :param conditions: Expression that must evaluate to True/False.

    :type conditions: numpy.ndarray

    :returns: numpy.ndarray -- Array of indexes where the condition is False.
    """

    return np.where(conditions == False)[0]
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
        print "         ", filename

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
            strarr = str.split(line)
            if len(strarr) > 0:
                tbl.append(strarr)

    return np.array(tbl, dtype='float64')
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def chunk(a, b, length=1, verbose=0):
    """
    Produces an array, of 2x1, delimiting ranges between 'a' and 'b', with
    the requested length.

    :param a: Minimum of the range to consider.

    :type a: float

    :param b: Maximum of the range to consider.

    :type b: float

    :param length: Length of the array to create. @CHASE - confirm description.@

    :type length: int @CHASE - Can this be a float?@

    :param verbose: @CHASE - This is not used, can it be removed?@

    :type verbose: int

    :returns: list -- The segmented array. @CHASE - please refine if needed.@
    """

    if not length:
        return [[a, b]]
    else:
        arr = np.array([np.arange(a, b, length), np.arange(a, b,
                                                           length)+length]).T
        if arr[-1][1] > b:
            arr[-1][1] = b
        return arr.tolist()
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def chunks(array, length=10, verbose=0):
    """
    Takes an array of ranges and produces a new array of ranges
    within the initial array that are of the requested length (or smaller).

    :param array: The array of ranges (2-element lists) to operate on.
    @CHASE - Confirm this description.@

    :type array: list @CHASE - Is this is list or numpy.ndarray?@

    :param length: Length of the arrays to create. @CHASE - confirm description@

    :type length: int @CHASE - int or float?@

    :param verbose: @CHASE - This parameter is not used in 'chunk' and can be
    removed.  If not, it should be a Boolean True/False and not an int.@

    :type verbose: int

    :returns: list -- The updated set of ranges. @CHASE - Is this a list or
    numpy.ndarray? Also please update the description.@
    """

    out = []

    for a in array:
        out += chunk(a[0], a[1], length=length, verbose=verbose)

    return out
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

    # @CHASE - Remove these two lines below?@
    # sep = np.sin(d0)*np.sin(d1)+np.cos(d0)*np.cos(d1)*np.cos((ra1-ra0)*dtor)
    # r = np.arccos(sep)*radeg

    a = np.sin((d2-d1)/2.)**2.+np.cos(d1)*np.cos(d2)*np.sin((r2-r1)/2.)**2.
    r = 2*np.arcsin(np.sqrt(a))

    # @CHASE - Remove these commented out lines below?@
    # zero = (np.isfinite(r) == False)
    # if any(zero):
    #	r[zero] = 0.0

    return r*ra2deg
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def intersect(r1, r2):
    """
    Returns the intersection of r1 and r2.
    @CHASE - There is a 'FIXME' here, is this method broken?@

    :param r1: The first set of values.

    :type r1: list @CHASE - confirm data type.@

    :param r2: The second set of values.

    :type r2: list @CHASE - confirm data types.@

    :returns: tuple -- A 2-element tuple containing the minimum and maximum
    value representing the intersection between the two sets of values.
    @CHASE - Confirm data type and description.@
    """

    # @CHASE - Here is the 'FIXME' line below.@
    # FIXME
    t0, t1 = np.array(r1)
    trange = np.array(r2)

    if (t0 <= trange[0]) & (t1 > trange[1]):
        return trange[0], trange[1]
    elif (t0 <= trange[0]) & (t1 < trange[1]):
        return trange[0], t1
    elif (t0 > trange[0]) & (t1 > trange[1]):
        return t0, trange[1]
    elif (t0 > trange[0]) & (t1 <= trange[1]):
        return t0, t1
    else:
        print t0, t1, trange
        return None
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def algebraicIntersection(steps, tranges):
    """
    Returns intervals defining the intersection of r1 and r2.

    :param steps: @CHASE - please describe this parameter.@

    :type steps: @CHASE - Is this is a list?@

    :param tranges: @CHASE - please describe this parameter.@

    :type tranges: @CHASE - Is this a list?@

    :returns: list -- The intervals defining the intersection.
    @CHASE - Please confirm data type and description.@
    """

    t0 = np.array(steps)[:, 0]
    t1 = np.array(steps)[:, 1]

    sect = []

    for trange in tranges:
        ix = np.where(
            ((t0 >= trange[0]) & (t0 < trange[1])) |
            ((t1 >= trange[0]) & (t1 < trange[1])) |
            ((t0 <= trange[0]) & (t1 > trange[1])))[0]
        sect += [intersect([t0[i], t1[i]], trange) for i in ix]

    return sect
# ------------------------------------------------------------------------------
