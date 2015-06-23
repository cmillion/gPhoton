import os as _os
from .. import cal_dir
from ..MCUtils import get_fits_data, get_fits_header, get_tbl_data

def check_band(band):
    if not band in ['NUV','FUV']:
        raise ValueError('Band must be NUV or FUV.')
    return band

def check_dimension(dimension):
    if not dimension in ['x','y']:
        raise ValueError('Dimension must be x or y.')
    return dimension

def read_data(fn):
    path = _os.path.join(cal_dir,fn)
    if '.fits' in fn:
        return get_fits_data(path), get_fits_header(path)
    elif '.tbl' in fn:
        return get_tbl_data(path)
    else:
        raise ValueError('Unrecognized data type: {ext}'.format(fn[-4:]))

def wiggle(band,dimension):
    fn = '{b}_wiggle_{d}.fits'.format(
            b=check_band(band),d=check_dimension(dimension))
    return read_data(fn)

def avgwalk(band,dimension):
    fn = '{b}_avgwalk_{d}.fits'.format(
            b=check_band(band),d=check_dimension(dimension))
    return read_data(fn)

def walk(band,dimension):
    fn = '{b}_walk_{d}.fits'.format(
            b=check_band(band),d=check_dimension(dimension))
    return read_data(fn)

def linearity(band,dimension):
    fn = '{b}_NLC_{d}_det2sky.fits'.format(
            b=check_band(band),d=check_dimension(dimension))
    return read_data(fn)

def flat(band):
    fn = '{b}_flat.fits'.format(b=check_band(band))
    return read_data(fn)

def distortion_filenames(band,dimension,eclipse,raw_stimsep):
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
            b=check_band(band).lower(),d=check_dimension(dimension),i=index)
    return read_data(fn)

def offset(band,dimension):
    fn = '{b}_d{d}_fdttdc_coef_0.tbl'.format(
            b=check_band(band).lower(),d=check_dimension(dimension))
    return read_data(fn)

def mask_filename(band,calpath):
    fn = '{b}_mask.fits'.format(b=check_band(band))
    return read_data(fn)
