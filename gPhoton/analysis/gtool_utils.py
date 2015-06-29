__version__ = '1.0'

"""
.. module:: gtool_utils

   :synopsis: A variety of utility functions, designed especially for 
   use with the gTool package.

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

from astropy.time import Time

#--------------------
def calculate_jd(galex_time):
    """
    Calculates the Julian date, in TDB time scale, given a GALEX time.

    :param galex_time: A GALEX timestamp.
    
    :type galex_time: float

    :returns: float -- The time converted to a Julian date, in the TDB 
    time scale.
    """
    
    """ Convert the GALEX timestamp to a Unix timestamp. """
    this_unix_time = Time(galex_time + 315964800., format="unix", 
                          scale="utc")

    """ Convert the Unix timestamp to a Julian date, measured in the 
    TDB scale. """
    this_jd_time = this_unix_time.tdb.jd

    return this_jd_time
#--------------------

#--------------------
def calculate_caldat(galex_time):
    """
    Calculates a Gregorian calendar date given a GALEX time, in the UTC 
    time scale.

    :param galex_time: A GALEX timestamp.
    
    :type galex_time: float

    :returns: float -- The time converted to a Gregorian calendar date, 
    in the UTC time scale.
    """
    
    """ Convert the GALEX timestamp to a Unix timestamp. """
    this_unix_time = Time(galex_time + 315964800., format="unix", 
                          scale="utc")

    """ Convert the Unix timestamp to a Julian date, measured in the 
    TDB scale. """
    this_caldat_time = this_unix_time.iso

    return this_caldat_time
#--------------------

