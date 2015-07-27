__version__ = '1.25.5'
import time
import os.path as _osp
pkg_dir = _osp.abspath(_osp.dirname(__file__))
cal_dir = _osp.join(pkg_dir, 'cal')
time_id = int(time.time()) # Define this session by time at import.
