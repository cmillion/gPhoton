from __future__ import absolute_import
# Core and Third Party imports.
import os.path as _osp
import time
__version__ = '1.28.3'
pkg_dir = _osp.abspath(_osp.dirname(__file__))
cal_dir = _osp.join(pkg_dir, 'cal')
time_id = int(time.time()) # Define this session by time at import.
# Import the main modules here to simplify syntax to just, e.g., gPhoton.gFind
# [DEVEL] NOTE: It is CRITICAL that these lines go AFTER the constants above.
#         Due to the (circular) import structure, setup.py fails unless these
#         come first, even though it's against PEP8 standard.
from gPhoton.gAperture import gaperture as gAperture
from gPhoton.gFind import gfind as gFind
from gPhoton.gMap import gmap as gMap
