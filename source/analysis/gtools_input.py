#!/usr/bin/env python

__version__ = '1.0'

"""
.. module:: gtools_input

   :synopsis: Given one or more coordinates and optional identifiers, 
   creates a list of gTarget class objects; the object class that is 
   used throughout gTools to send and store results.  Includes gTarget 
   read and write modules to save states to disk.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import argparse
import os
import pandas as pd

""" These module-level variables ensure both the ArgumentParser and 
module use the same defaults. """
#--------------------
ifile_default = None
coordtype_default = "radec"
galcoords_default = None
radeccoords_default = None
target_id_default = "Target_1"
output_dir_default = "StateFiles"
#--------------------

#--------------------
class gToolsInputError(Exception, object):
    """
    This class defines a generic Exception to use for errors raised in the gTools "input" module.  It simply prints the given string when raising the exception. e.g.,

    .. code-block:: python
     
         raise gToolsInputError("Print this string")
         gToolsInputError: *** gTools Input Error: 'Print this string'
    """

    def __init__(self, value):
        """
        Initiate the Exception.

        :param value: The string to print on error.

        :type value: str
        """
        self.value = value

    def __str__(self):
        """
        Overrides the str function for this class.
        """
        return "*** gTools Input Error: "+repr(self.value)
#--------------------

#--------------------
class gTarget(object):
    pass
#--------------------

#--------------------
def read_gTarget():
    pass
#--------------------

#--------------------
def write_gTarget():
    pass
#--------------------

#--------------------
def setup_args():
    """
    Set up command-line arguments and options.

    :returns: ArgumentParser -- Stores arguments and options.
    """
    parser = argparse.ArgumentParser(description="""Read in targets 
                                     for use with gTools.""")

    parser.add_argument("-f", action="store", type=str, dest="ifile", 
                        default=ifile_default, 
                        help="""Full path to input file containing 
                        coordinates (and optional identifiers) of 
                        targets to use.  Format of the file should be:
                        (ID), RA or GLON, DEC or GLAT, separated by 
                        commas.  Lines that begin with a pound sign (#) 
                        are considered comments and ignored.  Include 
                        the file name in the path.  If an input file 
                        is specified, the command-line coordinates and 
                        identifier are ignored.""")
    parser.add_argument("--coordtype", action="store", 
                        default=coordtype_default, 
                        choices=["radec", "galactic"], 
                        help="""What type of coordinates are supplied 
                        in the target list?  Default="%(default)s".""")

    coordinate_group = parser.add_mutually_exclusive_group()
    coordinate_group.add_argument("--galcoords", action="store", 
                                  nargs=2, default=galcoords_default, 
                                  help="""Galactic longitude and 
                                  lattitude of target, if not 
                                  supplying a list of targets.  
                                  Specify both in decimal degrees.  
                                  Ignored if a list of targets is 
                                  specified with the -f option.  
                                  Example: --galcoords 123.45 67.89""")
    coordinate_group.add_argument("--radeccoords", action="store", 
                                  nargs=2, default=radeccoords_default, 
                                  help="""RA and DEC of target, if not 
                                  supplying a list of targets.  
                                  Specify both in decimal degrees.  
                                  Ignored if a list of targets is 
                                  specified with the -f option.  
                                  Example: --radeccoords 123.45 67.89""")

    parser.add_argument("-i", action="store", dest="target_id", 
                        default=target_id_default, 
                        help="""[Optional] Identifier of target if not 
                        supplying a list of targets.  If no ID is 
                        provided, target will be referred to as 
                        "%(default)s".  Ignored if a list of targets 
                        is specified with the -f option.""")

    parser.add_argument("-o", action="store", dest="output_dir", 
                        default=output_dir_default, 
                        help="""[Optional] Full path where state files 
                        will be saved.  Default is to save in a folder
                         called "%(default)s" in the current working 
                        directory.""")

    return parser
#--------------------

#--------------------
def check_input_options(args):
    """
    Check that input arguments satisfy some minimum requirements.

    :param args: Stored arguments and options.

    :type args: argparse.Namespace object.

    :raises: gToolsInputError, IOError
    """

    """ Either an input file OR one of the coordinates must have been 
    specified. """
    if (args.ifile is None and args.galcoords is None and 
        args.radeccoords is None):
        raise gToolsInputError("Must specify either an input file or "
                               "set of coordinates.")

    """ If an input file is specified, ignore the other parameters by 
    making sure they are set to None. """
    if args.ifile is not None:
        args.galcoords = [None, None]
        args.radeccoords = [None, None]

    """ Make sure the input file exists. """
    if not os.path.isfile(args.ifile):
        raise IOError("File containing list of targets not found.  "
                      "Looking for " + args.ifile)
#--------------------

#--------------------
def input_targets(ifile=ifile_default, coordtype=coordtype_default, 
                  galcoords=galcoords_default, 
                  radeccoords=radeccoords_default, 
                  target_id=target_id_default, 
                  output_dir=output_dir_default):
    """
    Parses the input file or command line for target identifiers and 
    coordinates.  Creates gTarget object classes for each one, does
    initial population of information, and saves state files.
        
    :param ifile: Input file containing target coordinates and, 
    optionally, identifiers.

    :type ifile: str

    :param coordtype:  Does the input file contain galactic coordinates 
    (glon, glat) or equitorial coordinates (RA, DEC)?  Should be either 
    "galactic" or "radec".  Default is """ + coordtype_default + """.
    
    :type coordtype: str

    :param galcoords:  The galactic longitude and latitude for the 
    target, in the absence of an input list of targets.  Specify in 
    decimal degrees.

    :type galcoords: list

    :param radeccoords:  The RA and DEC for the target, in the absence 
    of an input list of targets.  Specify in degrees.

    :type radeccoords: list

    :param target_id:  Optional identifier for the target, in the 
    absence of an input list of targets.  If not supplied, a default 
    value of """ + taregt_id_default + """ will be used.

    :type target_id: str

    :param output_dir:  Optional output path for state files.  By 
    default, state files are stored in a subdirectory 
    called """ + output_dir_default + """ in the current working 
    directory.

    :type output_dir: str

    :returns: list -- A list of gTarget objects for each target.

    :raises: gToolsInputError
    """

    """ Read in list of targets, or, create a one-row data frame from
    the input coordinates. """
    if ifile is not None:
        file_contents = pd.io.parsers.read_csv(ifile, skipinitialspace=
                                               True, header=None, 
                                               comment='#', 
                                               skip_blank_lines=True)
        n_cols = len(file_contents.columns)
        n_rows = len(file_contents.index)
        if n_cols != 2 and n_cols !=3:
            raise gToolsInputError('Input file must have exactly 2 or '
                                   '3 columns: (ID), RA/GLON, DEC/GLAT.'
                                   '  Found ' + str(n_cols) + 
                                   ' colunns.')

        if n_cols == 2:
            """ Then we need to add our own target IDs. """
            file_contents.insert(0, "ID", ["Target_{0:02d}".format(x+1)
                                           for x in xrange(n_rows)])
        """ Now there are three columns.  Make sure the headers are set
        in order. """
        file_contents.rename(columns={file_contents.columns[0]:"ID", 
                                      file_contents.columns[1]:"COORD1",
                                      file_contents.columns[2]:"COORD2"
                                      }, inplace=True)
    elif radeccoords is not None:
        """ Create one-row data frame, don't convert coordinats yet. """
        pass
    elif galcoords is not None:
        """ Create one-row data frame, don't convert coordinats yet. """
        pass
#--------------------

#--------------------
if __name__ == "__main__":
    """ Create ArgumentParser object that holds arguments and 
    options. """
    args = setup_args().parse_args()

    """ Check arguments and options. """
    check_input_options(args)

    """ Call primary method. """
    input_targets(ifile=args.ifile, coordtype=args.coordtype, 
                  galcoords=args.galcoords, 
                  radeccoords=args.radeccoords, 
                  target_id=args.target_id, output_dir=args.output_dir)
#--------------------
