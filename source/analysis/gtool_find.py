__version__ = '1.0'

"""
.. module:: gtool_find

   :synopsis: Given a gTarget state file, will execute gFind to 
   determine the exposure time information for the given target.  Serves
   as a wrapper to gFind itself.  Also does some additional calculations
   to convert reported GALEX times into JD and Gregorian dates and times
   .

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import argparse
import os
from gtool_read import read as gt_read

""" <DEVEL> Note that this hack to make it so that the user can import `stitch_components` directly as a module or run it from the command line as __main__ has the side effect of importing this module twice, despite my best efforts to work around it.  I don't think it will be a major issue, but worth thinking about in the future. </DEVEL> """
if __package__ is None:
    import sys, os
    gtfind_dir = os.path.dirname(os.path.abspath(__file__))
    gtfind_pardir = os.path.dirname(gtfind_dir)
    sys.path.insert(1, gtfind_pardir)
    from gFind import gFind
#    __package__ = str("source")
#    __name__ = str(__package__+"."+__name__)

""" These module-level variables ensure both the ArgumentParser and 
module use the same defaults. """
#--------------------
ifile_default = None
#--------------------

#--------------------
def setup_args():
    """
    Set up command-line arguments and options.

    :returns: ArgumentParser -- Stores arguments and options.
    """
    parser = argparse.ArgumentParser(description="""Read in state files 
                                     for use with gtool_find.""")

    parser.add_argument("ifile", action="store", type=str, 
                        default=ifile_default, 
                        help="""Full path and file name of the input 
                        state file.""")
    return parser
#--------------------

#--------------------
def check_input_options(args):
    """
    Check that input arguments satisfy some minimum requirements.

    :param args: Stored arguments and options.

    :type args: argparse.Namespace object.

    :raises: IOError
    """

    """ Make sure the input file exists. """
    if args.ifile is not None:
        if not os.path.isfile(args.ifile):
            raise IOError("State file not found.  Looking for " + 
                          args.ifile)
#--------------------

#--------------------
def gtool_find(ifile=ifile_default):
    """
    Reads in the state file and executes a gFind command to obtain
    exposure time information.  Updates the state file with this 
    information.
    
    :param ifile: Input state file, it will also be the file that is
    updated with the results of the gFind query.

    :type ifile: str

    """

    """ Read in the state file. """
    state_file_content = gt_read(ifile)

    """ Test call of gFind. """
    exptime = gFind(band="both", exponly=True, skypos=[state_file_content.ra, state_file_content.dec])
    print exptime["NUV"]["expt"]
#--------------------

#--------------------
if __name__ == "__main__":
    """ Create ArgumentParser object that holds arguments and 
    options. """
    args = setup_args().parse_args()

    """ Check arguments and options. """
    check_input_options(args)

    """ Call primary method. """
    gtool_find(args.ifile)
#--------------------


