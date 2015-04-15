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
def json_encoder(obj):
    """
    Defines a method to use when serializing into JSON.
    """
    return obj.__dict__
#--------------------

#--------------------
def write(gt, odir, prefix):
    """
    Writes the gTarget (or list of gTarget) object(s) to a JSON file in
    the specified output directory.

    :param gt: A single gTarget object, or a list of gTarget objects, 
    to write to output file(s).

    :type gt: gTarget object or list

    :param odir: Output directory where the files should be written.  
    If the directory does not exist, it will be created.

    :type odir: str

    :param prefix: Prefix string to prepend to output state files.

    :type prefix: str

    :returns: list -- A list of the output files written to disk.
    """

    """ Create the output directory if it does not yet exists. """
    if not os.path.isdir(odir):
        os.mkdir(odir)

    """ Make sure output directory includes the trailing slash. """
    odir = os.path.join(odir,'')

    """ If the input object is not a list, put it into one. """
    if not isinstance(gt, list):
        gt = [gt]

    """ Create a list of output files. """
    output_files = [odir+prefix+"_"+t.id+".json" for t in gt]

    """ For each object in the list, output to a JSON file. """
    for t,ofile in zip(gt,output_files):
        with open(ofile, 'wb') as of:
            json.dump(t,of,default=json_encoder,indent=4)

    """ Return the list of output files. """
    return output_files
#--------------------

