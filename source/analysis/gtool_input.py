__version__ = '1.0'

"""
.. module:: gtool_input

   :synopsis: Given one or more coordinates and optional identifiers, 
   creates a list of gTarget class objects; the object class that is 
   used throughout gTool to send and store results.  Includes gTarget 
   read and write modules to save states to disk.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import argparse
import os
import numpy as np
import pandas as pd
from astropy import units
from astropy.coordinates import SkyCoord
from gtool_write import write as gt_write

""" These module-level variables ensure both the ArgumentParser and 
module use the same defaults. """
#--------------------
ifile_default = None
coordtype_default = None
galcoords_default = None
radeccoords_default = None
target_id_default = "Target_1"
output_dir_default = "StateFiles"
output_prefix_default = "gTool"
#--------------------

#--------------------
class gToolInputError(Exception, object):
    """
    This class defines a generic Exception to use for errors raised in 
    the gTool "input" module.  It simply prints the given string when 
    raising the exception. e.g.,

    .. code-block:: python
     
         raise gToolInputError("Print this string")
         gToolInputError: *** gTool Input Error: 'Print this string'
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
        return "*** gTool Input Error: "+repr(self.value)
#--------------------

#--------------------
class gTarget(object):
    """
    This class defines a single target within gPhoton/gTool.  It keeps
    track of the target identifier, coordinates, locations of gPhoton
    output files, location of diagnostic plots, lightcurve parameters,
    aperture sizes, etc.
    """

    def __init__(self, targetid, ra, dec, glon, glat, 
                 fuv_tot_exptime=None, nuv_tot_exptime=None,
                 fuv_timerange_start=None, fuv_timerange_end=None,
                 nuv_timerange_start=None, nuv_timerange_end=None,
                 fuv_start_stop=None, nuv_start_stop=None):
        """
        :param targetid: ID/name of the target.

        :type targetid: str
        
        :param ra: Right Ascension of the target in degrees.
        
        :type ra: numpy.float64

        :param dec: Declination of the target in degrees.
        
        :type dec: numpy.float64

        :param glon: Galactic longitude of the target in degrees.

        :type glon: numpy.float64

        :param glat: Galactic lattitude of the target in degrees.
        
        :type glat: numpy.float64
        """
        self.id = targetid
        self.ra = ra
        self.dec = dec
        self.glon = glon
        self.glat = glat
        self.fuv_tot_exptime = fuv_tot_exptime
        self.nuv_tot_exptime = nuv_tot_exptime
        self.fuv_timerange_start = fuv_timerange_start
        self.fuv_timerange_end = fuv_timerange_end
        self.nuv_timerange_start = nuv_timerange_start
        self.nuv_timerange_end = nuv_timerange_end
        self.fuv_start_stop = fuv_start_stop
        self.nuv_start_stop = nuv_start_stop

    @classmethod
    def from_json(self, idict):
        """
        Create a gTarget object from a dictionary of values.

        :param idict: Dictionary containing values to populate the 
        gTarget object with.

        :type idict: dict
        """
        
        """ Make sure the dict contains all the expected keywords. """
        if (set([u'id', u'ra', u'dec', u'glon', u'glat', 
                 u'fuv_tot_exptime', u'nuv_tot_exptime', 
                 u'fuv_timerange_start', u'nuv_timerange_start', 
                 u'fuv_timerange_end', u'nuv_timerange_end', 
                 u'fuv_start_stop', u'nuv_start_stop']) <= 
            set(idict.keys())):
            return gTarget(idict[u"id"], idict[u"ra"], idict[u"dec"], 
                           idict[u"glon"], idict[u"glat"], 
                           fuv_tot_exptime = idict[u"fuv_tot_exptime"], 
                           nuv_tot_exptime = idict[u"nuv_tot_exptime"], 
                           fuv_timerange_start = 
                           idict[u"fuv_timerange_start"], 
                           fuv_timerange_end = 
                           idict[u"fuv_timerange_end"], 
                           nuv_timerange_start = 
                           idict[u"nuv_timerange_start"],
                           nuv_timerange_end = 
                           idict[u"nuv_timerange_end"], 
                           fuv_start_stop = idict[u"fuv_start_stop"], 
                           nuv_start_stop = idict[u"nuv_start_stop"])
        else:
            raise gToolInputError("Input dict does not have expected "
                                  "set of keys.")            
#--------------------

#--------------------
def setup_args():
    """
    Set up command-line arguments and options.

    :returns: ArgumentParser -- Stores arguments and options.
    """
    parser = argparse.ArgumentParser(description="""Read in targets 
                                     for use with gTool.""")

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
                        identifier are ignored.  RA/DEC coordinates are 
                        assumed to be in the ICRS frame.""")
    parser.add_argument("--coordtype", action="store", 
                        default=coordtype_default, 
                        choices=["radec", "galactic"], 
                        help="""What type of coordinates are supplied 
                        in the target list?  RA/DEC are assumed to be 
                        in the ICRS frame.""")

    coordinate_group = parser.add_mutually_exclusive_group()
    coordinate_group.add_argument("--galcoords", action="store", 
                                  nargs=2, default=galcoords_default, 
                                  help="""Galactic longitude and 
                                  lattitude of the target, if not 
                                  supplying a list of targets.  
                                  Specify both in decimal degrees.  
                                  Ignored if a list of targets 
                                  is specified with the -f option.  
                                  Example: --galcoords 123.45 67.89""")
    coordinate_group.add_argument("--radeccoords", action="store", 
                                  nargs=2, default=radeccoords_default, 
                                  help="""RA and DEC of the target, if 
                                  not supplying a list of targets.  
                                  Specify both in decimal degrees.  
                                  The ICRS frame is used.  
                                  Ignored if a list of targets is 
                                  specified with the -f option.  
                                  Example: --radeccoords 123.45 67.89"""
                                  )

    parser.add_argument("-i", action="store", dest="target_id", 
                        default=target_id_default, 
                        help="""[Optional] Identifier of target if not 
                        supplying a list of targets.  If no ID is 
                        provided, target will be referred to as 
                        "%(default)s".  Ignored if a list of targets 
                        is specified with the -f option.""")

    parser.add_argument("-n", action="store", dest="output_prefix", 
                        default=output_prefix_default, 
                        help="""[Optional] Prefix to use for the output
                        state file.  Default is "%(default)s".""")

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

    :raises: gToolInputError, IOError
    """

    """ Either an input file OR one of the coordinates must have been 
    specified. """
    if (args.ifile is None and args.galcoords is None and 
        args.radeccoords is None):
        raise gToolInputError("Must specify either an input file or "
                               "set of coordinates.")

    """ If an input file is specified, ignore the other parameters by 
    making sure they are set to None. """
    if args.ifile is not None:
        args.galcoords = [None, None]
        args.radeccoords = [None, None]

    """ Make sure the input file exists. """
    if args.ifile is not None:
        if not os.path.isfile(args.ifile):
            raise IOError("File containing list of targets not found.  "
                          "Looking for " + args.ifile)
#--------------------

#--------------------
def input_targets(ifile=ifile_default, coordtype=coordtype_default, 
                  galcoords=galcoords_default, 
                  radeccoords=radeccoords_default, 
                  target_id=target_id_default, 
                  output_dir=output_dir_default,
                  output_prefix=output_prefix_default):
    """
    Parses the input file or command line for target identifiers and 
    coordinates.  Creates gTarget object classes for each one, does
    initial population of information, and saves state files.
        
    :param ifile: Input file containing target coordinates and, 
    optionally, identifiers.

    :type ifile: str

    :param coordtype:  Does the input file contain galactic coordinates 
    (glon, glat) or equitorial coordinates (RA, DEC)?  Should be either 
    "galactic" or "radec".
    
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
    value of """ + target_id_default + """ will be used.

    :type target_id: str

    :param output_dir:  Optional output path for state files.  By 
    default, state files are stored in a subdirectory 
    called """ + output_dir_default + """ in the current working 
    directory.

    :type output_dir: str

    :param output_prefix: Prefix to use on the output state files.  
    Default is """ + output_prefix_default + """.

    :type output_prefix: str

    :returns: list -- A list of gTarget objects for each target.

    :raises: gToolInputError
    """

    """ Read in list of targets, or, create a one-row data frame from
    the input coordinates. """
    if ifile is not None:
        """ The user must specify a coordinate type. """
        if coordtype is None:
            raise gToolInputError('You must specify a coordinate '
                                   'type when supplying a target list.')

        file_contents = pd.io.parsers.read_csv(ifile, skipinitialspace=
                                               True, header=None, 
                                               comment='#', 
                                               skip_blank_lines=True)
        n_cols = len(file_contents.columns)
        n_rows = len(file_contents.index)
        if n_cols != 2 and n_cols !=3:
            raise gToolInputError('Input file must have exactly 2 or '
                                   '3 columns: (ID), RA/GLON, DEC/GLAT.'
                                   '  Found ' + str(n_cols) + 
                                   ' colunns.')

        if n_cols == 2:
            """ Then we need to add our own target IDs. """
            file_contents.insert(0, "ID", [("Target_{0:0" +
                                            str(len(str(n_rows))) + 
                                            "d}").format(x+1)
                                           for x in xrange(n_rows)])
        """ Now there are three columns.  Make sure the headers are set
        in order. """
        file_contents.rename(columns={file_contents.columns[0]:"ID", 
                                      file_contents.columns[1]:"COORD1",
                                      file_contents.columns[2]:"COORD2"
                                      }, inplace=True)
    elif radeccoords is not None:
        """ Make the coordtype "radec". """
        coordtype = "radec"
        """ Create one-row data frame, don't convert coordinats yet. """
        file_contents = pd.DataFrame([{"ID":target_id,
                                       "COORD1":radeccoords[0],
                                       "COORD2":radeccoords[1]}])
        n_rows = 1
        n_cols = 3
    elif galcoords is not None:
        """ Make the coordtype "galactic". """
        coordtype = "galactic"
        """ Create one-row data frame, don't convert coordinats yet. """
        file_contents = pd.DataFrame([{"ID":target_id,
                                       "COORD1":galcoords[0],
                                       "COORD2":galcoords[1]}])
        n_rows = 1
        n_cols = 3

    """ Make sure the COORD columns are 64-bit floats and not strings 
    (they can be strings if given from the command line. """
    file_contents["COORD1"] = file_contents["COORD1"].astype('float64')
    file_contents["COORD2"] = file_contents["COORD2"].astype('float64')

    """ Now assign the Coordinate columns in the DataFrame to the 
    appropriate column name (RA, GLON, etc.), and add the other set of 
    coordinates. """
    if coordtype == "radec":
        file_contents.rename(columns={"COORD1":"RA",
                                      "COORD2":"DEC"}, inplace=True)
        """ Add the galactic coordinates. """
        gal_coords = SkyCoord(ra=file_contents["RA"]*units.degree,
                              dec=file_contents["DEC"]*units.degree,
                              frame="icrs").galactic
        file_contents.insert(len(file_contents.columns), "GLON", 
                             [x.l.deg for x in gal_coords])
        file_contents.insert(len(file_contents.columns), "GLAT", 
                             [x.b.deg for x in gal_coords])
    elif coordtype == "galactic":
        file_contents.rename(columns={"COORD1":"GLON",
                                      "COORD2":"GLAT"}, inplace=True)
        radec_coords = SkyCoord(l=file_contents["GLON"]*units.degree,
                              b=file_contents["GLAT"]*units.degree,
                              frame="galactic").icrs
        file_contents.insert(len(file_contents.columns), "RA", 
                             [x.ra.deg for x in radec_coords])
        file_contents.insert(len(file_contents.columns), "DEC", 
                             [x.dec.deg for x in radec_coords])

    """ Create a list of gTarget objects. """
    gtargets_list = [gTarget(file_contents["ID"][i],
                             file_contents["RA"][i],
                             file_contents["DEC"][i],
                             file_contents["GLON"][i],
                             file_contents["GLAT"][i]
                             ) for i in xrange(n_rows)]

    """ Write each gTarget object to disk. """
    return gt_write(gtargets_list, output_dir, output_prefix)
#--------------------

#--------------------
if __name__ == "__main__":
    """ Create ArgumentParser object that holds arguments and 
    options. """
    args = setup_args().parse_args()

    """ Check arguments and options. """
    check_input_options(args)

    """ Call primary method. """
    output_files = input_targets(ifile=args.ifile, 
                                 coordtype=args.coordtype, 
                                 galcoords=args.galcoords, 
                                 radeccoords=args.radeccoords, 
                                 target_id=args.target_id, 
                                 output_dir=args.output_dir,
                                 output_prefix=args.output_prefix)
#--------------------
