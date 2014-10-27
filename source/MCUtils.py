# A library of generic utilities programs that Chase wants to keep
# separate in order to re-use across projects.
from sys import stdout
import numpy as np
import requests
import time
from astropy.io import fits as pyfits
import csv

def error(program,note=''):
	"""Print ERROR with callback and note."""
        print '*** ERROR: '+str(program)
        if note:
                print '    '+str(note)
        return

def area(radius):
    """Returns area of a circle of given radius."""
    return np.pi*radius**2.

def distance(a,b,c,d):
    """Computes Euclidean distance between [a,b] and [c,d]."""
    return np.sqrt( (a-c)**2. + (b-d)**2. )

def rotvec(vector,theta):
	"""Rotate vectors clockwise by theta degrees
	vector -- must have shape (2,n)
	theta  -- must have shape (n,)
	returns w/ shape (2,n)
	"""
	R = np.array([[np.cos(np.radians(theta)),-np.sin(np.radians(theta))],
		      [np.sin(np.radians(theta)), np.cos(np.radians(theta))]])
	return np.dot(R.T,vector)[0]

def rms(data):
        return np.sqrt(np.square(np.array(data)-np.array(data).mean()).mean())

def print_inline(text,blanks=60):
	"""For updating text in place without a carriage return."""
        stdout.write(" "*blanks+"\r")
        stdout.write(str(str(text)+'\r'))
        stdout.flush()
        return

def manage_requests(query,maxcnt=100,wait=10):
	""" Make simple 'requests' calls slightly more robust against network
    issues.
    """
	cnt = 0
	while cnt < maxcnt:
		try:
			r = requests.get(query)
			# HACK: This specifically tests for the return structure
			#  from the MAST photon database and is therefore
			#  not a good general test condition.
			test = r.json()['Tables'][0]['Rows']
			return r
		except:
			time.sleep(wait)
			cnt += 1
			print_inline('Query retry attempt '+str(int(cnt))+'.')


	print 'Query unsuccessful after '+str(int(maxcnt))+' attempts.'
	print '		'+str(query)
	return False

def wheretrue(conditions):
	"""Returns indices for which the input conditions are true.
	Conditions must be a numpy array.
	Example:
	> a = np.array([1,2,3,4,5,6,7,8,9])
	> wheretrue(a>5)
	array([5, 6, 7, 8])
	> a[wheretrue(a>5)]
	array([6, 7, 8, 9])
	"""
	return np.where(conditions==True)[0]

# Returns an array of indexes for which the input conditions are false.
# - as above, conditions must be a numpy array.
# Example: as above
def wherefalse(conditions):
	"""Returns indices for which the input conditions are false.
	See the docstring for wheretrue() for an example.
	"""
	return np.where(conditions==False)[0]

def find_nearest_lower(array,value):
	"""Finds the index of the value in the array that is closest
    without going over.

	This function assumes that:
		1. value is within the range of array
		2. array is ordered
		3. array has no gaps
	"""
        idx = (np.abs(array-value)).argmin()
        if array[idx]>value:
                idx -= 1
        return idx

def get_fits_data(filename,dim=0,verbose=0):
	"""Reads FITS data. A wrapper for common pyfits commands."""
	if verbose:
        	print "         ",filename
        hdulist = pyfits.open(filename,memmap=1)
        data = hdulist[dim].data
        hdulist.close()
        return data

def get_fits_header(filename):
	"""Reads a FITS header. A wrapper for common pyfits commands."""
        hdulist = pyfits.open(filename,memmap=1)
        htab = hdulist[0].header
        hdulist.close()
        return htab

def get_tbl_data(filename,comment='|'):
	"""Reads data from a table into a numpy array."""
        f = open(filename)
        lines = f.readlines()
        tbl = []
        for line in lines:
                if line[0] != comment:
                        strarr = str.split(line)
                        if len(strarr) > 0:
                                tbl.append(strarr)

        return np.array(tbl,dtype='float64')

def chunk(a,b,length=1,verbose=0):
	"""Produces an array of 2x1 delimiting ranges between a and b with
    the requested length.
    """
	if not length:
		return [[a,b]]
	else:
		arr=np.array([np.arange(a,b,length),np.arange(a,b,length)+length]).T
		if arr[-1][1]>b:
			arr[-1][1]=b
		return arr.tolist()

def chunks(array,length=10,verbose=0):
	"""Takes an array of ranges and produces a new array of ranges
    within the initial array that are of the requested length (or smaller).
    """
	out = []
	for a in array:
		out += chunk(a[0],a[1],length=length,verbose=verbose)
	return out

def angularSeparation(ra1,dec1,ra2,dec2):
	"""Compute angular separation in degrees of points on the sky.
    It is important, especially for small angular separations, that the
    values for ra[01] and dec[01] have precision of float64 or better.
	Update: Now uses the haversine formula which is stable for small angles.
    """
	d2r = np.pi/180.
	ra2deg = 1./d2r
	d1 = dec1*d2r
	d2 = dec2*d2r
	r1 = ra1*d2r
	r2 = ra2*d2r
	#sep = np.sin(d0)*np.sin(d1)+np.cos(d0)*np.cos(d1)*np.cos((ra1-ra0)*dtor)
	#r = np.arccos(sep)*radeg
	a = np.sin((d2-d1)/2.)**2.+np.cos(d1)*np.cos(d2)*np.sin((r2-r1)/2.)**2.
	r = 2*np.arcsin(np.sqrt(a))
	#zero = (np.isfinite(r) == False)
	#if any(zero):
	#	r[zero] = 0.0
	return r*ra2deg


def intersect(r1,r2):
    #FIXME
    t0,t1=np.array(r1)
    trange=np.array(r2)
    """Returns the intersection of r1 and r2."""
    if (t0<=trange[0]) & (t1>trange[1]):
        return trange[0],trange[1]
    elif (t0<=trange[0]) & (t1<trange[1]):
        return trange[0],t1
    elif (t0>trange[0]) & (t1>trange[1]):
        return t0,trange[1]
    elif (t0>trange[0]) & (t1<=trange[1]):
        return t0,t1
    else:
        print t0,t1,trange
        return None

def algebraicIntersection(steps,tranges):
    """Returns intervals defining the intersection of r1 and r2."""
    t0=np.array(steps)[:,0]
    t1=np.array(steps)[:,1]
    sect = []
    for trange in tranges:
        ix = np.where(((t0>=trange[0]) & (t0<trange[1])) | ((t1>=trange[0]) & (t1<trange[1])) | ((t0<=trange[0]) & (t1>trange[1])))[0]
        sect+=[intersect([t0[i],t1[i]],trange) for i in ix]
    return sect
