__version__ = '1.0'

"""
.. module:: gtool_write

   :synopsis: Writes a gTarget object to a file on disk.  This 
   enables a user to save their analysis of a target at any stage, i.e.,
   it serves as a "state" file.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import json
import os

#--------------------
def json_default(obj):
    """
    Defines a default method to use when serializing into JSON.
    """
    return obj.__dict__
#--------------------

#--------------------
def write(gt, odir):
    """
    Writes the gTarget (or list of gTarget) object(s) to a JSON file in
    the specified output directory.

    :param gt: A single gTarget object, or a list of gTarget objects, 
    to write to output file(s).

    :type gt: gTarget object or list

    :param odir: Output directory where the files should be written.  
    If the directory does not exist, it will be created.

    :type odir: str
    """

    """ Create the output directory if it does not yet exists. """
    if not os.path.isdir(odir):
        os.mkdir(odir)

    """ Make sure output directory includes the trailing slash. """
    odir = os.path.join(odir,'')

    """ If the input object is not a list, put it into one. """
    if not isinstance(gt, list):
        gt = [gt]

    """ For each object in the list, output to a JSON file. """
    for t in gt:
        with open(odir+"gTool_"+t.id+".json", 'wb') as of:
            json.dump(t,of,default=json_default,indent=4)
#--------------------

