# This file defines common queries that are passed to the GALEX photon
# database at MAST.
from MCUtils import manage_requests
import CalUtils
from galextools import isPostCSP
from gPhoton import time_id

tscale = 1000. # time in the database is "integer-ized" by multiplying by this

"""The following three global variables are used in constructing a properly
formatted query to the MAST database. Don't change them unless you know what
you're doing.
"""
baseURL = 'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query='
baseDB = 'GPFCore.dbo'
MCATDB = 'GR6Plus7.dbo'
formatURL = ' -- '+str(time_id)+'&format=json&timeout={}'

def hasNaN(query):
    """Check if there is NaN in a query (or any string) and, if so, raise an
    exception because that probably indicates that something has gone wrong.
    """
    if 'NaN' in query:
        raise RuntimeError("Malformed query: contains NaN values.")
    return

def getValue(query,verbose=0,retries=20):
    """Manage a database call which returns a single value."""
    hasNaN(query)
    if verbose>2:
        print query
    try:
        out = float(manage_requests(query,
                            maxcnt=retries).json()['Tables'][0]['Rows'][0][0])
    except:
        raise# RuntimeError('Connection timeout.')
    return out

def getArray(query,verbose=0,retries=20):
    """Manage a database call which returns an array of values."""
    hasNaN(query)
    if verbose>2:
        print query
    try:
        out = manage_requests(query,maxcnt=retries).json()['Tables'][0]['Rows']
    except:
        raise# RuntimeError('Connection timeout.')
    return out

def mcat_sources(band,ra0,dec0,radius,maglimit=20):
    ''' Return the MCAT _coadd_ sources given sky position and search radius
    (and optional lower magnitude limit).
    Columns are:
    [0,RA],[1,Dec],[2,NUV_mag],[3,FUV_mag],[4,FoV_radius],[5,NUV_skybg],
    [6,FUV_skybg],[7,NUV_FWHM_world],[8,FUV_FWHM_world],
    [9:15,FUV_mag_aper_1:7],[16:22,NUV_mag_aper_1:7]
    [23:29,FUV_magerr_aper_1:7],[30:36,NUV_magerr_aper1:7]
    '''
    # 1=nuv, 2=fuv, 3=both
    bandflag = 1 if band=='NUV' else 2
    # fGetNearbyObjEq takes radius in arcminutes
    # TODO: Add exposure time.
    return (str(baseURL)+
        'select ra, dec, nuv_mag, fuv_mag, fov_radius, nuv_skybg, fuv_skybg,'
        ' nuv_fwhm_world, fuv_fwhm_world, fuv_mag_aper_1, fuv_mag_aper_2,'
        ' fuv_mag_aper_3, fuv_mag_aper_4, fuv_mag_aper_5, fuv_mag_aper_6,'
        ' fuv_mag_aper_7, nuv_mag_aper_1, nuv_mag_aper_2, nuv_mag_aper_3,'
        ' nuv_mag_aper_4, nuv_mag_aper_5, nuv_mag_aper_6, nuv_mag_aper_7'
        ' fuv_magerr_aper_1, fuv_magerr_aper_2, fuv_magerr_aper_3,'
        ' fuv_magerr_aper_4, fuv_magerr_aper_5, fuv_magerr_aper_6,'
        ' fuv_magerr_aper_7, nuv_magerr_aper_1, nuv_magerr_aper_2,'
        ' nuv_magerr_aper_3, nuv_magerr_aper_4, nuv_magerr_aper_5,'
        ' nuv_magerr_aper_6, nuv_magerr_aper_7'
        ' from '+str(MCATDB)+'.photoobjall as p inner join '+str(MCATDB)+
        '.photoextract as pe on p.photoextractid=pe.photoextractid inner join '+
        str(MCATDB)+'.fgetnearbyobjeq('+repr(float(ra0))+', '+
        repr(float(dec0))+', '+
        str(radius*60.)+') as nb on p.objid=nb.objid and (band=3 or band='+
        str(bandflag)+') and '+str(band)+'_mag<'+str(maglimit)+
        str(formatURL))

def mcat_visit_sources(ra0,dec0,radius):
    ''' Return the MCAT per-visit sources given sky position and search radius.
    The columns are as follows:
    [0,objid],[1,ra],[2,dec],[3,NUV_mag],[4,FUV_mag],[5,FoV_radius],
    [6,NUV_skybg],[7,FUV_skybg],[8,NUV_FWHM],[9,FUV_FWHM],[10,FUV_expt],
    [11,NUV_expt],[12:18,FUV_mag_aper_1:7],[19:25,NUV_mag_aper_1:7],
    [26:32,FUV_magerr_aper_1:7],[33:39,NUV_magerr_aper_1:7]
    '''
    return (str(baseURL)+
        'select vpo.objid, ra, dec, nuv_mag, fuv_mag, fov_radius, nuv_skybg,'
        ' fuv_skybg, nuv_fwhm_world, fuv_fwhm_world, vpe.fexptime,'
        ' vpe.nexptime, fuv_mag_aper_1, fuv_mag_aper_2, fuv_mag_aper_3,'
        ' fuv_mag_aper_4, fuv_mag_aper_5, fuv_mag_aper_6, fuv_mag_aper_7,'
        ' nuv_mag_aper_1, nuv_mag_aper_2, nuv_mag_aper_3, nuv_mag_aper_4,'
        ' nuv_mag_aper_5, nuv_mag_aper_6, nuv_mag_aper_7,'
        ' fuv_magerr_aper_1, fuv_magerr_aper_2, fuv_magerr_aper_3,'
        ' fuv_magerr_aper_4, fuv_magerr_aper_5, fuv_magerr_aper_6,'
        ' fuv_magerr_aper_7, nuv_magerr_aper_1, nuv_magerr_aper_2,'
        ' nuv_magerr_aper_3, nuv_magerr_aper_4, nuv_magerr_aper_5,'
        ' nuv_magerr_aper_6, nuv_magerr_aper_7'
        ' from '+str(MCATDB)+'.visitphotoobjall as vpo'
        ' inner join '+str(MCATDB)+'.visitphotoextract'
        ' as vpe on vpo.photoextractid=vpe.photoextractid inner join'
        ' '+str(MCATDB)+'.fGetNearbyVisitObjEq('+repr(float(ra0))+','+
        repr(float(dec0))+', '+str(radius*60.)+
        ') as nb on vpo.objid=nb.objid'+str(formatURL))

def mcat_objid_search(objid,mode='visit'):
    """Return a bunch of observation data for a visit level objid (ggoid).
    Doing the same for coadd level data is not yet supported.
    """
    return (str(baseURL)+
        'select objid, minPhotoObsDate, maxPhotoObsDate, obs_date, obsdatim,'
        ' nobssecs, fobssecs, nexptime, fexptime, nexpstar, nexpend, fexpstar,'
        ' fexpend from '+str(MCATDB)+'.visitphotoobjall as vp inner join '+
        str(MCATDB)+'.imgrun as ir on vp.photoextractid=ir.imgrunid'
        ' inner join '+str(MCATDB)+'.visitphotoextract as vpe on'
        ' vp.photoextractid=vpe.photoextractid where objid = '+
        str(long(objid))+str(formatURL))

def exposure_ranges(band,ra0,dec0,t0=1,t1=10000000000000,detsize=1.25):
    """Returns a list of times (in one second increments) where data exists
    with an aspect solution within detsize of [ra0,dec0].
    """
    # If band is set to False (or equivalent), search on both bands
    if not band:
        band = 'FUV/NUV'
    return (str(baseURL)+
        'select distinct time from '+str(baseDB)+
        '.fGetNearbyAspectEq('+repr(float(ra0))+','+repr(float(dec0))+',(('+
            str(detsize)+'/2.0)*60.0),'+
            str(long(t0*tscale))+','+str(long(t1*tscale))+')'
        ' where band=\''+str(band)+'\' or band=\'FUV/NUV\' order by time'+
        str(formatURL))

def exposure_range(band,ra0,dec0,t0=1,t1=10000000000000):
    """Find time ranges for which data exists at a given position."""
    return (str(baseURL)+
        'select startTimeRange, endTimeRange'
        ' from '+str(baseDB)+
            '.fGetTimeRanges('+str(int(t0))+','+str(int(t1))+','+
            repr(float(ra0))+','+repr(float(dec0))+') where band=\''+
            str(band)+'\''+str(formatURL))

def aperture(band,ra0,dec0,t0,t1,radius):
    """Integrate counts over an aperture at a position."""
    return (str(baseURL)+
        'Select  sum(photonCount)'
        ' from '+str(baseDB)+'.fGetNearbyObjEqCount'+str(band)+
        '('+repr(float(ra0))+','+repr(float(dec0))+
        ','+str(radius)+','+str(long(t0*tscale))+','+
        str(long(t1*tscale))+',0)'+str(formatURL))

def deadtime1(band,t0,t1):
    """Return the global counts of non-NULL data."""
    return ('{baseURL}select count(*) from {baseDB}.{band}PhotonsV where '+
            'time >= {t0} and time < {t1}{formatURL}').format(baseURL=baseURL,
                baseDB=baseDB, band=band, t0=str(long(t0*tscale)),
                t1=str(long(t1*tscale)), formatURL=formatURL)

def deadtime2(band,t0,t1):
    """Return the global counts for NULL data."""
    return ('{baseURL}select count(*) from {baseDB}.{band}PhotonsNULLV where '+
            'time >= {t0} and time < {t1}{formatURL}').format(baseURL=baseURL,
                baseDB=baseDB, band=band, t0=str(long(t0*tscale)),
                t1=str(long(t1*tscale)), formatURL=formatURL)

def deadtime(band,t0,t1,feeclkratio=0.966,tec2fdead=5.52e-6):
    """Return the emperically determined deadtime correction based upon the
    global count rate.
    """
    scale = tec2fdead/feeclkratio
    return (str(baseURL)+
        'select sum(dt)*'+'%0.30f'%scale+' / ('+repr(t1)+'-'+repr(t0)+')'
        ' from(select count(*) as dt from '+str(baseDB)+'.'+str(band)+
        'PhotonsNULLV where time >= '+str(long(t0*tscale))+' and time < '+
        str(long(t1*tscale))+' union all select count(*) as dt from '+
        str(baseDB)+'.'+str(band)+
        'PhotonsV where time >= '+str(long(t0*tscale))+' and time < '+
        str(long(t1*tscale))+') x'+str(formatURL))

def globalcounts(band,t0,t1):
    """Return the total number of detector events within a time range."""
    return ('{baseURL}select count(t) from (select time as t from {baseDB}.{band}PhotonsV '+
            'where time >= {t0} and time < {t1} union all select time as t '+
            'from {baseDB}.{band}PhotonsNULLV where time >= '+
            '{t0} and time < {t1}) x{formatURL}').format(baseURL=baseURL,
                baseDB=baseDB, band=band, t0=str(long(t0*tscale)),
                t1=str(long(t1*tscale)), formatURL=formatURL)


def alltimes(band,t0,t1):
    """Return the time stamps of every detector event within a time range."""
    return ('{baseURL}select t from (select time as t from {baseDB}.{band}PhotonsV '+
            'where time >= {t0} and time < {t1} union all select time as t '+
            'from {baseDB}.{band}PhotonsNULLV where time >= '+
            '{t0} and time < {t1}) x{formatURL}').format(baseURL=baseURL,
                baseDB=baseDB, band=band, t0=str(long(t0*tscale)),
                t1=str(long(t1*tscale)), formatURL=formatURL)

def uniquetimes(band,t0,t1):
    """Return the _unique_ timestamps for all detector events within range."""
    #return ('{baseURL}select distinct t from (select time as t from #{baseDB}.{band}PhotonsV '+
#            'where time >= {t0} and time < {t1} union all select time as t '+
#            'from {baseDB}.{band}PhotonsNULLV where time >= '+
#            '{t0} and time < {t1}) x order by t'+
#            '{formatURL}').format(baseURL=baseURL,
#                baseDB=baseDB, band=band, t0=str(long(t0*tscale)),
#                t1=str(long(t1*tscale)), formatURL=formatURL)
    return ('{baseURL}select distinct time from {baseDB}.{band}PhotonsV '+
            'where time >= {t0} and time < {t1} order by time'+
            '{formatURL}').format(baseURL=baseURL,
                baseDB=baseDB, band=band, t0=str(long(t0*tscale)),
                t1=str(long(t1*tscale)), formatURL=formatURL)

def boxcount(band,t0,t1,xr,yr):
    """Find the number of events inside of a box defined by [xy] range in
    detector space coordinates. This is useful for pulling out stim events.
    """
    return (str(baseURL)+'select count(*) from '+str(baseDB)+'.'+str(band)+
        'PhotonsNULLV where time >= '+str(long(t0*tscale))+' and time < '+
        str(long(t1*tscale))+' and x >= '+str(xr[0])+' and x < '+str(xr[1])+
        ' and y >= '+str(yr[0])+' and y < '+str(yr[1])+str(formatURL))

def detbox(band,t0,t1,xr,yr):
    """Return all events inside a box defined in detector space by [xy]
    range. Created as a sanity check for stim events."""
    return ('{baseURL}select time,x,y,ya from {baseDB}.{band}PhotonsNULLV '+
            'where time >= {t0} and time < {t1} and '+
            'x >= {xmin} and x < {xmax} and y >= {ymin} and y < {ymax}'+
            '{formatURL}').format(baseURL=baseURL, baseDB=baseDB, band=band,
                t0=str(long(t0*tscale)), t1=str(long(t1*tscale)),
                xmin=xr[0], xmax=xr[1], ymin=yr[0], ymax=yr[1],
                formatURL=formatURL)

# FIXME: convert t0 to eclipse
def stimcount(band,t0,t1,margin=[90.01,90.01],aspum=68.754932/1000.,
              eclipse=None):
    """Return stim counts."""
    if not eclipse:
        eclipse = 55000 if isPostCSP(t0) else 30000
    if isPostCSP(t0):
        margin[1]=180.02
    avgstim = CalUtils.avg_stimpos(band,eclipse)
    return ('{baseURL}select count(*) from {baseDB}.{band}PhotonsNULLV '+
            'where time >= {t0} and time < {t1} and ('+
            '((x >= {x10} and x < {x11}) and (y >= {y10} and y < {y11})) or '+
            '((x >= {x20} and x < {x21}) and (y >= {y20} and y < {y21})) or '+
            '((x >= {x30} and x < {x31}) and (y >= {y30} and y < {y31})) or '+
            '((x >= {x40} and x < {x41}) and (y >= {y40} and y < {y41}))'+
            '){formatURL}').format(baseURL=baseURL, baseDB=baseDB, band=band,
                t0=str(long(t0*tscale)), t1=str(long(t1*tscale)),
                x10=(avgstim['x1']-margin[0])/aspum,
                x11=(avgstim['x1']+margin[0])/aspum,
                y10=(avgstim['y1']-margin[1])/aspum,
                y11=(avgstim['y1']+margin[1])/aspum,
                x20=(avgstim['x2']-margin[0])/aspum,
                x21=(avgstim['x2']+margin[0])/aspum,
                y20=(avgstim['y2']-margin[1])/aspum,
                y21=(avgstim['y2']+margin[1])/aspum,
                x30=(avgstim['x3']-margin[0])/aspum,
                x31=(avgstim['x3']+margin[0])/aspum,
                y30=(avgstim['y3']-margin[1])/aspum,
                y31=(avgstim['y3']+margin[1])/aspum,
                x40=(avgstim['x4']-margin[0])/aspum,
                x41=(avgstim['x4']+margin[0])/aspum,
                y40=(avgstim['y4']-margin[1])/aspum,
                y41=(avgstim['y4']+margin[1])/aspum,formatURL=formatURL)


def boxcentroid(band,t0,t1,xr,yr):
    """Find the mean position of events inside of a box in detector space."""
    return (str(baseURL)+
        'select avg(x), avg(y) from '+str(baseDB)+
        '.'+str(band)+'PhotonsNULLV where time >= '+
        str(long(t0*tscale))+' and time < '+str(long(t1*tscale))+
        ' and x >= '+str(xr[0])+' and x < '+str(xr[1])+' and y >= '+
        str(yr[0])+' and y < '+str(yr[1])+str(formatURL))

def boxtimes(band,t0,t1,xr,yr):
    """Get the list of times for events inside of a box in detector space."""
    return (str(baseURL)+
        'select time from '+str(baseDB)+'.'+str(band)+
        'PhotonsNULLV where time >= '+
        str(long(t0*tscale))+' and time < '+str(long(t1*tscale))+
        ' and x >= '+
        str(xr[0])+' and x < '+str(xr[1])+' and y >= '+str(yr[0])+' and y < '+
        str(yr[1])+str(formatURL))

def centroid(band,ra0,dec0,t0,t1,radius):
    return (str(baseURL)+
        'select avg(ra), avg(dec) from '+str(baseDB)+'.'+str(band)+
        'PhotonsV where time >= '+
        str(long(t0*tscale))+' and time < '+str(long(t1*tscale))+
        ' and ra >= '+
        repr(ra0-radius)+' and ra < '+repr(ra0+radius)+' and dec >= '+
        repr(dec0-radius)+' and dec < '+repr(dec0+radius)+str(formatURL))

def allphotons(band,ra0,dec0,t0,t1,radius):
    """Grab the major columns for all events within an aperture."""
    return (str(baseURL)+
        'select time,ra,dec,xi,eta,x,y from '+str(baseDB)+'.fGetNearbyObjEq'+
        str(band)+'AllColumns('+repr(float(ra0))+','+repr(float(dec0))+','+
        repr(radius)+','+
        str(long(t0*tscale))+','+str(long(t1*tscale))+',0)'+str(formatURL))

# Shutter correction
#  i.e. number of 0.05s gaps in data
def shutter(band,t0,t1):
    return (str(baseURL)+
        'select shutter*0.05 from '+str(baseDB)+'.fGet'+str(band)+'Shutter('+
        str(long(t0*tscale))+','+str(long(t1*tscale))+')'+str(formatURL))

def shutdead(band,t0,t1):
    tt0, tt1 = [long(t*tscale) for t in [t0, t1]]
    return ('{baseURL}SELECT shutter*0.05 FROM {baseDB}'
        '.fGet{band}Shutter({tt0},{tt1}) AS '
        'time UNION ALL SELECT SUM(dt) * 0.0000057142857142857145 / '
        '({t1}-{t0}) AS dead FROM(SELECT count(*) AS dt FROM {baseDB}'
        '.{band}PhotonsNULLV WHERE time >= {tt0} AND time < {tt1} UNION ALL '
        'SELECT count(*) AS dt FROM {baseDB}.{band}PhotonsV WHERE time '
        '>= {tt0} AND time < {tt1}) x{fmt}'.format(baseURL=baseURL,
                band=band.upper(), tt0=tt0, tt1=tt1, t0=t0, t1=t1,
                baseDB=baseDB, fmt=formatURL))


def exptime(band,t0,t1,stepsz=1.):
    return (str(baseURL)+
        'select * from '+str(baseDB)+'.fGet'+str(band)+'EffectiveExposureTime('+
        str(long(t0*tscale))+','+str(long(t1*tscale))+','+str(stepsz)+')'+
        str(formatURL))

def aspect(t0,t1):
    """Return the aspect information based on time range."""
    return (str(baseURL)+
        'select eclipse, filename, time, ra, dec, twist, flag, ra0, dec0,'
        ' twist0 from aspect where time >= '+str(long(t0*tscale))+
        ' and time < '+str(long(t1*tscale))+' order by time'+str(formatURL))

# Return the aspect information based on eclipse
def aspect_ecl(eclipse):
    """Return the aspect information based upon an eclipse number."""
    return (str(baseURL)+
        'select eclipse, filename, time, ra, dec, twist, flag, ra0, dec0,'
        ' twist0 from aspect where eclipse='+str(eclipse)+' order by time'+
        str(formatURL))

def aspect_skypos(ra,dec,detsize=1.25):
    """Return the aspect information based upon sky position and det radius."""
    return (str(baseURL)+
        "select eclipse, filename, time, ra, dec, twist, flag, ra0, dec0,"
        " twist0 from aspect where ra >= "+repr(ra-detsize/2.)+
        " and ra < "+repr(ra+detsize/2.)+" and dec >= "+repr(dec-detsize/2.)+
        " and dec < "+repr(dec+detsize/2.)+" order by time"+str(formatURL))

# Return data within a box centered on ra0, dec0 with sides of length 2*radius
# TODO: deprecate this and rename it skybox()
def box(band,ra0,dec0,t0,t1,radius):
    return (str(baseURL)+
        'select time,ra,dec from '+str(baseDB)+'.'+str(band)+
        'PhotonsV where time >= '+
        str(long(t0*tscale))+' and time < '+str(long(t1*tscale))+
        ' and ra >= '+
        repr(ra0-radius)+' and ra < '+repr(ra0+radius)+' and dec >= '+
        repr(dec0-radius)+' and dec < '+repr(dec0+radius)+' and flag=0'+
        str(formatURL))

# Return data within a rectangle centered on ra0, dec0
# TODO: deprecate this and rename it skyrect()
def rect(band,ra0,dec0,t0,t1,ra,dec):
    #time,ra,dec,xi,eta,x,y
    return (str(baseURL)+
        'select time,ra,dec,xi,eta,x,y from '+str(baseDB)+'.fGetObjFromRect'+str(band)+'('+
        repr(ra0-ra/2.)+','+repr(ra0+ra/2.)+','+repr(dec0-dec/2.)+','+
        repr(dec0+dec/2.)+','+str(long(t0*tscale))+','+
        str(long(t1*tscale))+',0)'+str(formatURL))

def skyrect(band,ra0,dec0,t0,t1,ra,dec):
    return '{baseURL}select time,ra,dec,xi,eta,x,y from {baseDB}.fGetObjFromRect{band}AllColumns({ra_min},{ra_max},{dec_min},{dec_max},{t0},{t1},0){formatURL}'.format(baseURL=baseURL,baseDB=baseDB,band=band,ra_min=repr(ra0-ra/2.),ra_max=repr(ra0+ra/2.),dec_min=repr(dec0-dec/2.),dec_max=repr(dec0+dec/2.),t0=str(long(t0*tscale)),t1=str(long(t1*tscale)),formatURL=formatURL)
