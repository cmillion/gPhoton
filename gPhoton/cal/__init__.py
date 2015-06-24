import os as _os
from .. import cal_dir
from ..MCUtils import get_fits_data, get_fits_header, get_tbl_data

def check_band(band):
    if not band in ['NUV','FUV']:
        raise ValueError('Band must be NUV or FUV')
    return band

def check_xy(xy):
    if not xy in ['x','y']:
        raise ValueError('xy must be x or y.')
    return xy

def read_data(fn,dim=0):
    path = _os.path.join(cal_dir,fn)
    if '.fits' in fn:
        return get_fits_data(path,dim=dim), get_fits_header(path)
    elif '.tbl' in fn:
        return get_tbl_data(path)
    else:
        raise ValueError('Unrecognized data type: {ext}'.format(fn[-4:]))

def wiggle(band,xy):
    fn = '{b}_wiggle_{d}.fits'.format(b=check_band(band),d=check_xy(xy))
    return read_data(fn)

def wiggle2():
    """The post-CSP wiggle file."""
    return read_data('WIG2_Sep2010.fits',dim=1)

def avgwalk(band,xy):
    fn = '{b}_avgwalk_{d}.fits'.format(
            b=check_band(band),d=check_xy(xy))
    return read_data(fn)

def walk(band,xy):
    fn = '{b}_walk_{d}.fits'.format(
            b=check_band(band),d=check_xy(xy))
    return read_data(fn)

def walk2():
    """The post-CSP walk file."""
    return read_data('WLK2_Sep2010.fits',dim=1)

def clock2():
    """The post-CSP clock file."""
    return read_data('CLK2_Sep2010.fits',dim=1)

def linearity(band,xy):
    fn = '{b}_NLC_{d}_det2sky.fits'.format(
            b=check_band(band),d=check_xy(xy))
    return read_data(fn)

def flat(band):
    fn = '{b}_flat.fits'.format(b=check_band(band))
    return read_data(fn)

def distortion(band,xy,eclipse,raw_stimsep):
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
            b=check_band(band).lower(),d=check_xy(xy),i=index)
    return read_data(fn)

def offset(band,xy):
    fn = '{b}_d{d}_fdttdc_coef_0.tbl'.format(
            b=check_band(band).lower(),d=check_xy(xy))
    return read_data(fn)

def mask(band):
    fn = '{b}_mask.fits'.format(b=check_band(band))
    return read_data(fn)
