__version__ = '1.0'

"""
.. module:: gtool_read

   :synopsis: Reads in a gTarget object from a file on disk.  This 
   enables a user to restart analysis of a target at any stage, i.e.,
   it serves as a "state" file.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import json
import os
from gtool_input import gTarget

#--------------------
def json_decoder(obj):
    """
    Defines a method to use when deflating from JSON.
    """
    return gTarget.from_json(obj)
#--------------------

#--------------------
def read(ifile):
    """
    Reads in the gTarget information contained in the JSON input file 
    and returns the object.

    :param ifile: The file name of the JSON object file to read 
    (including the full path).

    :type ifile: str

    :returns: gTarget -- A gTarget object generated from the information
    contained within the JSON file.
    """

    """ Check that the file exists. """
    if not os.path.isfile(ifile):
        raise IOError("Input file not found.  Looking for " + ifile + 
                      ".")

    """ Read in the file's contents. """
    with open(ifile, 'rb') as input_file:
        this_gTarget = json.load(input_file, 
                                         object_hook=json_decoder)

    return this_gTarget
#--------------------

#--------------------
if __name__ == "__main__":
    """
    If called from the command line, simply reads in the gTarget JSON 
    file and returns True/False if the returned object is of type 
    gTarget.
    """
    import sys
    if len(sys.argv) != 2:
        raise ValueError("Unexpected number of arguments encountered."
                         "  Usage: python gtool_read.py <file>")
    else:
        ifile = sys.argv[1]

    gt = read_gTarget(ifile)

    print "Returned object is gTarget object: " + str(isinstance(gt, 
                                                    gTarget))
#--------------------
