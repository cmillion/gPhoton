from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
from builtins import range
import sys
import os as _os
import numpy as np
from astropy.io import fits

# gPhoton imports.
from gPhoton.MCUtils import get_fits_data, get_fits_header, get_tbl_data
from gPhoton import cal_dir
#from gPhoton.download import download_with_progress_bar

# Remote repository for GALEX calibration files.
cal_url = 'https://archive.stsci.edu/prepds/gphoton/cal/cal/'

"""
The following three functions are substantially derived from code in
https://github.com/astroML/astroML and so carry the following license:

Copyright (c) 2012-2013, Jacob Vanderplas All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

if (sys.version_info[0] == 3):
    from urllib.request import urlopen
    from urllib.error import HTTPError
    from urllib.parse import urlencode
    from io import BytesIO
else:
    from urllib2 import urlopen
    from urllib2 import HTTPError
    from urllib import urlencode
    from cStringIO import StringIO as BytesIO

def url_content_length(fhandle):
    if (sys.version_info[0] == 3):
        length = dict(fhandle.info())['Content-Length']
    else:
        length = fhandle.info().getheader('Content-Length')
    return int(length.strip())

def bytes_to_string(nbytes):
    if nbytes < 1024:
        return '%ib' % nbytes
    nbytes /= 1024.
    if nbytes < 1024:
        return '%.1fkb' % nbytes
    nbytes /= 1024.
    if nbytes < 1024:
        return '%.2fMb' % nbytes
    nbytes /= 1024.
    return '%.1fGb' % nbytes

def download_with_progress_bar(data_url, file_path):
    if not _os.path.exists(_os.path.dirname(file_path)):
        _os.makedirs(_os.path.dirname(file_path))
    num_units = 40
    fhandle = urlopen(data_url)
    content_length = url_content_length(fhandle)
    chunk_size = content_length // num_units
    print('Downloading {url} to {cal_dir}'.format(
                            url=data_url,cal_dir=_os.path.dirname(file_path)))
    nchunks = 0
    buf = BytesIO()
    content_length_str = bytes_to_string(content_length)
    while True:
        next_chunk = fhandle.read(chunk_size)
        nchunks += 1
        if next_chunk:
            buf.write(next_chunk)
            s = ('[' + nchunks * '='
                + (num_units - 1 - nchunks) * ' '
                + ']  %s / %s   \r' % (bytes_to_string(buf.tell()),
                                        content_length_str))
        else:
            sys.stdout.write('\n')
            break

        sys.stdout.write(s)
        sys.stdout.flush()

    buf.seek(0)
    open(file_path, 'wb').write(buf.getvalue())
    return

### END code derived from astroML

def check_band(band):
    if not band in ['NUV', 'FUV']:
        raise ValueError('Band must be NUV or FUV')
    return band

def check_xy(xy):
    if not xy in ['x', 'y']:
        raise ValueError('xy must be x or y.')
    return xy

def read_data(fn, dim=0):
    path = _os.path.join(cal_dir, fn)
    # Download the file if it doesn't exist locally.
    if not _os.path.exists(path):
        data_url='{b}/{f}'.format(b=cal_url,f=fn)
        fitsdata = download_with_progress_bar(data_url, path)
    if '.fits' in fn:
        return get_fits_data(path, dim=dim), get_fits_header(path)
    elif '.tbl' in fn:
        return get_tbl_data(path)
    else:
        raise ValueError('Unrecognized data type: {ext}'.format(ext=fn[-4:]))

def wiggle(band, xy):
    fn = '{b}_wiggle_{d}.fits'.format(b=check_band(band), d=check_xy(xy))
    return read_data(fn)

def wiggle2():
    """The post-CSP wiggle file."""
    return read_data('WIG2_Sep2010.fits', dim=1)

def avgwalk(band, xy):
    fn = '{b}_avgwalk_{d}.fits'.format(
        b=check_band(band), d=check_xy(xy))
    return read_data(fn)

def walk(band, xy):
    fn = '{b}_walk_{d}.fits'.format(
        b=check_band(band), d=check_xy(xy))
    return read_data(fn)

def walk2():
    """The post-CSP walk file."""
    return read_data('WLK2_Sep2010.fits', dim=1)

def clock2():
    """The post-CSP clock file."""
    return read_data('CLK2_Sep2010.fits', dim=1)

def linearity(band, xy):
    fn = '{b}_NLC_{d}_det2sky.fits'.format(
        b=check_band(band), d=check_xy(xy))
    return read_data(fn)

def addbuffer(fn):
# Adds a 1 pixel buffer around all masked (==0) regions of a map.
    m, h = read_data(fn)
    ix = np.where(m == 0)
    for i in range(-1, 2):
        for j in range(-1, 2):
            try:
                m[ix[0]+i, ix[1]+j] = 0
            except IndexError:
                continue
    return m, h

def flat(band, buffer=False):
    fn = '{b}_flat.fits'.format(b=check_band(band))
    return addbuffer(fn) if buffer else read_data(fn)

def distortion(band, xy, eclipse, raw_stimsep):
    index = ''
    if band == 'NUV':
        if eclipse > 37460:
            if raw_stimsep < 5136.3:
                index = 'a'
            elif raw_stimsep < 5137.25:
                index = 'b'
            else:
                index = 'c'
    fn = '{b}_distortion_cube_d{d}{i}.fits'.format(
        b=check_band(band).lower(), d=check_xy(xy), i=index)
    return read_data(fn)

def offset(xy):
    fn = 'fuv_d{d}_fdttdc_coef_0.tbl'.format(d=check_xy(xy))
    return read_data(fn)

def mask(band, buffer=False):
    fn = '{b}_mask.fits'.format(b=check_band(band))
    return addbuffer(fn) if buffer else read_data(fn)
