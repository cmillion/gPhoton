# A library of generic utilities programs that Chase wants to keep
# separate in order to re-use across projects.
from sys import stdout
import numpy as np
import requests
import time
from astropy.io import fits as pyfits
import csv

# Print an error with callback and note
def error(program,note=''):
        print '*** ERROR: '+str(program)
        if note:
                print '    '+str(note)
        return

def area(radius):
	return np.pi*radius**2.

# Rotate vectors clockwise by theta degrees
#  vector must have shape (2,n) and thata must have shape (n,)
#  will return shape (2,n)
def rotvec(vector,theta):
	R = np.array([[np.cos(np.radians(theta)),-np.sin(np.radians(theta))],
		      [np.sin(np.radians(theta)), np.cos(np.radians(theta))]])
	return np.dot(R.T,vector)[0]

def rms(data):
        return np.sqrt(np.square(np.array(data)-np.array(data).mean()).mean())

# Used for updating text output in place without doing a carriage return.
def print_inline(text):
        stdout.write("                                                                       \r")
        stdout.write(str(text+'\r'))
        stdout.flush()
        return

# A tool to make simple requests calls more robust against network issues.
def manage_requests(query,maxcnt=20,wait=1):
	cnt = 0
	while cnt < maxcnt:
		try:
			r = requests.get(query)
			# This specifically tests for the return structure
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
	return 0

# A tool to make simple grequest calls more robust against network issues.
def manage_grequests(queries,maxcnt=20,wait=1):
	cnt = 0
	while cnt < maxcnt:
		try:
			r = (grequests.get(q) for q in queries)
			return r
		except:
			time.sleep(wait)
			cnt +=1
			print_inline('Queries retry attempt'+str(int(cnt))+'.')

	print 'Queries unsuccessful after '+str(int(maxcnt))+' attempts.'
	for q in queries:
		print '		'+str(q)
	return 0

# Returns an array of indexes for which the input conditions are true.
# - conditions must be a numpy array.
# Example:
# >>> a = np.array([1,2,3,4,5,6,7,8,9])
# >>> wheretrue(a>5)
# array([5, 6, 7, 8])
# >>> a[wheretrue(a>5)]
# array([6, 7, 8, 9])
def wheretrue(conditions):
	return np.where(conditions==True)[0]

# Returns an array of indexes for which the input conditions are false.
# - as above, conditions must be a numpy array.
# Example: as above
def wherefalse(conditions):
	return np.where(conditions==False)[0]

def find_nearest_lower(array,value):
        # This assumes that:
        #   1. value is within the range of array
        #   2. array is ordered
        #   3. array has no gaps
        idx = (np.abs(array-value)).argmin()
        if array[idx]>value:
                idx -= 1
        return idx

def get_fits_data(filename,dim=0,verbose=1):
	if verbose:
        	print "         ",filename
        hdulist = pyfits.open(filename,memmap=1)
        data = hdulist[dim].data
        hdulist.close()
        return data

def get_fits_header(filename):
        hdulist = pyfits.open(filename,memmap=1)
        htab = hdulist[0].header
        hdulist.close()
        return htab

def get_tbl_data(filename,comment='|'):
        f = open(filename)
        lines = f.readlines()
        tbl = []
        for line in lines:
                if line[0] != comment:
                        strarr = str.split(line)
                        if len(strarr) > 0:
                                tbl.append(strarr)

        return np.array(tbl,dtype='float64')

# Produces an array of 2x1 delimiting ranges between a and b
#  with the requested length
def chunk(a,b,length=1,verbose=0):
	if not length:
		return [[a,b]]
	else:
		arr=np.array([np.arange(a,b,length),np.arange(a,b,length)+length]).T
		if arr[-1][1]>b:
			arr[-1][1]=b
		return arr.tolist()

def chunks(array,length=10,verbose=0):
	out = []
	for a in array:
		out += chunk(a[0],a[1],length=length,verbose=verbose)
	return out

# Generalized function to read csv files into python dictionary structures
#  with keys == column header names. Pass an array of column titles to
#  'columns'. Any entries of None or False will not be read.
# TODO: Add a way to specify numpy data types per column
# It seems that pandas can do this way better, so this function will go away soon.
def readcsv(csvfile,columns=None,verbose=0):
	data = {}
	# If the column names are not given, figure them out.
	# TODO: Actually pull out column names from the first line
	if not columns:
		# Name all the columns after their number
		columns = map(lambda x: str(x),
			range(len(csv.reader(open(csvfile,'rb')).next())))
	for i,col in enumerate(columns):
		if verbose>2:
			print i, col
		if col: 
			# TODO: Deal with commented lines
			data[col]=map(lambda row: row[i],
				      csv.reader(open(csvfile,'rb')))

	return data
