from __future__ import absolute_import
import time
import os.path as _osp
# Import the main modules here to simplify syntax to just, e.g., gPhoton.gFind
from gFind import gfind as gFind
from gAperture import gaperture as gAperture
from gMap import gmap as gMap

__version__ = '1.27.2'
pkg_dir = _osp.abspath(_osp.dirname(__file__))
cal_dir = _osp.join(pkg_dir, 'cal')
time_id = int(time.time()) # Define this session by time at import.
