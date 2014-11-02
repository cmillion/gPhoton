# This file defines common queries that are passed to the GALEX photon
# database at MAST.
from MCUtils import manage_requests
import CalUtils

"""The following three global variables are used in constructing a properly
formatted query to the MAST database. Don't change them unless you know what
you're doing.
"""
baseURL = 'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query='
baseDbo = 'Gr6plus7.dbo'
formatURL = '&format=json&timeout={}'

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
        raise RuntimeError('Connection timeout.')
    return out

def getArray(query,verbose=0,retries=20):
    """Manage a database call which returns an array of values."""
    hasNaN(query)
    if verbose>2:
        print query
    try:
        out = manage_requests(query,maxcnt=retries).json()['Tables'][0]['Rows']
    except:
        raise RuntimeError('Connection timeout.')
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
        ' from '+str(baseDbo)+'.photoobjall as p inner join '+str(baseDbo)+
        '.photoextract as pe on p.photoextractid=pe.photoextractid inner join '+
        str(baseDbo)+'.fgetnearbyobjeq('+str(ra0)+', '+str(dec0)+', '+
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
        ' from '+str(baseDbo)+'.visitphotoobjall as vpo'
        ' inner join '+str(baseDbo)+'.visitphotoextract'
        ' as vpe on vpo.photoextractid=vpe.photoextractid inner join'
        ' '+str(baseDbo)+'.fGetNearbyVisitObjEq('+str(ra0)+', '+str(dec0)+
        ', '+str(radius*60.)+') as nb on vpo.objid=nb.objid'+
        str(formatURL))

def mcat_objid_search(objid,mode='visit'):
    """Return a bunch of observation data for a visit level objid (ggoid).
    Doing the same for coadd level data is not yet supported.
    """
    return (str(baseURL)+
        'select objid, minPhotoObsDate, maxPhotoObsDate, obs_date, obsdatim,'
        ' nobssecs, fobssecs, nexptime, fexptime, nexpstar, nexpend, fexpstar,'
        ' fexpend from '+str(baseDbo)+'.visitphotoobjall as vp inner join '+
        str(baseDbo)+'.imgrun as ir on vp.photoextractid=ir.imgrunid'
        ' inner join '+str(baseDbo)+'.visitphotoextract as vpe on'
        ' vp.photoextractid=vpe.photoextractid where objid = '+
        str(long(objid))+str(formatURL))

def exposure_ranges(band,ra0,dec0,t0=1,t1=10000000000000,tscale=1000.,detsize=1.25):
    """Returns a list of times (in one second increments) where data exists
    with an aspect solution within detsize of [ra0,dec0].
    """
    # If band is set to False (or equivalent), search on both bands
    if not band:
        band = 'FUV/NUV'
    return (str(baseURL)+
        'select distinct time from'
        ' fGetNearbyAspectEq('+str(ra0)+','+str(dec0)+',(('+
            str(detsize)+'/2.0)*60.0),'+
            str(long(t0*tscale))+','+str(long(t1*tscale))+')'
        ' where band=\''+str(band)+'\' or band=\'FUV/NUV\' order by time'+
        str(formatURL))

def exposure_range(band,ra0,dec0,t0=1,t1=10000000000000):
    """Find time ranges for which data exists at a given position."""
    return (str(baseURL)+
        'select startTimeRange, endTimeRange'
        ' from fGetTimeRanges('+str(int(t0))+','+str(int(t1))+','+
            repr(ra0)+','+repr(dec0)+') where band=\''+str(band)+'\''+
        str(formatURL))

def aperture(band,ra0,dec0,t0,t1,radius,tscale=1000.):
    """Integrate counts over an aperture at a position."""
    return (str(baseURL)+
        'Select  sum(photonCount)'
        ' from dbo.fGetNearbyObjEqCount'+str(band)+'('+repr(ra0)+','+repr(dec0)+
            ','+str(radius)+','+str(long(t0*tscale))+','+
            str(long(t1*tscale))+',0)'+
        str(formatURL))

def deadtime1(band,t0,t1,tscale=1000.):
    """Return the global counts of non-NULL data."""
    return (str(baseURL)+
        'select count(*) from '+str(band)+'PhotonsV where time between '+
        str(long(t0*tscale))+' and '+str(long(t1*tscale))+str(formatURL))

def deadtime2(band,t0,t1,tscale=1000.):
    """Return the global counts for NULL data."""
    return (str(baseURL)+
        'select count(*) from '+str(band)+'PhotonsNULLV where time between '+
        str(long(t0*tscale))+' and '+str(long(t1*tscale))+str(formatURL))

def deadtime(band,t0,t1,tscale=1000.):
    """Return the emperically determined deadtime correction based upon the
    global count rate.
    """
    return (str(baseURL)+
        'select sum(dt)*0.0000057142857142857145 / ('+repr(t1)+'-'+repr(t0)+')'
        ' from(select count(*) as dt from '+str(band)+
        'PhotonsNULLV where time between '+str(long(t0*tscale))+' and '+
        str(long(t1*tscale))+' union all select count(*) as dt from '+
        str(band)+'PhotonsV where time between '+str(long(t0*tscale))+' and '+
        str(long(t1*tscale))+') x'+str(formatURL))

def boxcount(band,t0,t1,xr,yr,tscale=1000.):
    """Find the number of events inside of a box defined by [xy] range in
    detector space coordinates. This is useful for pulling out stim events.
    """
    return (str(baseURL)+'select count(*) from '+str(band)+
        'PhotonsNULLV where time between '+str(long(t0*tscale))+' and '+
        str(long(t1*tscale))+' and x between '+str(xr[0])+' and '+str(xr[1])+
        ' and y between '+str(yr[0])+' and '+str(yr[1])+str(formatURL))

# FIXME: convert t0 to eclipse
def stimcount(band,t0,t1,margin=90.01,aspum=68.754932/1000.,tscale=1000,eclipse=31000):
    """Return stim counts."""
    avgstim = CalUtils.avg_stimpos(band,eclipse)
    return (str(baseURL)+
        'select count(*) from '+str(band)+'PhotonsNULLV where time between '+
        str(long(t0*tscale))+' and '+str(long(t1*tscale))+' and ((x between '+
        str((avgstim['x1']-margin)/aspum)+' and '+
        str((avgstim['x1']+margin)/aspum)+' and y between '+
        str((avgstim['y1']-margin)/aspum)+' and '+
        str((avgstim['y1']+margin)/aspum)+') or (x between '+
        str((avgstim['x2']-margin)/aspum)+' and '+
        str((avgstim['x2']+margin)/aspum)+' and y between '+
        str((avgstim['y2']-margin)/aspum)+' and '+
        str((avgstim['y2']+margin)/aspum)+') or (x between '+
        str((avgstim['x3']-margin)/aspum)+' and '+
        str((avgstim['x3']+margin)/aspum)+' and y between '+
        str((avgstim['y3']-margin)/aspum)+' and '+
        str((avgstim['y3']+margin)/aspum)+') or (x between '+
        str((avgstim['x4']-margin)/aspum)+' and '+
        str((avgstim['x4']+margin)/aspum)+' and y between '+
        str((avgstim['y4']-margin)/aspum)+' and '+
        str((avgstim['y4']+margin)/aspum)+'))'+
        str(formatURL))

def boxcentroid(band,t0,t1,xr,yr,tscale=1000.):
    """Find the mean position of events inside of a box in detector space."""
    return (str(baseURL)+
        'select avg(x), avg(y) from NUVPhotonsNULLV where time between '+
        str(long(t0*tscale))+' and '+str(long(t1*tscale))+
        ' and x between '+str(xr[0])+' and '+str(xr[1])+' and y between '+
        str(yr[0])+' and '+str(yr[1])+str(formatURL))

def boxtimes(band,t0,t1,xr,yr,tscale=1000.):
    """Get the list of times for events inside of a box in detector space."""
    return (str(baseURL)+
        'select time from NUVPhotonsNULLV where time between '+
        str(long(t0*tscale))+' and '+str(long(t1*tscale))+' and x between '+
        str(xr[0])+' and '+str(xr[1])+' and y between '+str(yr[0])+' and '+
        str(yr[1])+str(formatURL))

def centroid(band,ra0,dec0,t0,t1,radius,tscale=1000.):
    return (str(baseURL)+
        'select avg(ra), avg(dec) from NUVPhotonsV where time between '+
        str(long(t0*tscale))+' and '+str(long(t1*tscale))+' and ra between '+
        repr(ra0-radius)+' and '+repr(ra0+radius)+' and dec between '+
        repr(dec0-radius)+' and '+repr(dec0+radius)+str(formatURL))

def allphotons(band,ra0,dec0,t0,t1,radius,tscale=1000.):
    """Grab the major columns for all events within an aperture."""
    return (str(baseURL)+
        'select time,ra,dec,xi,eta from dbo.fGetNearbyObjEq'+str(band)+
        'AllColumns('+repr(ra0)+','+repr(dec0)+','+repr(radius)+','+
        str(long(t0*tscale))+','+str(long(t1*tscale))+',0)'+str(formatURL))

# Shutter correction
#  i.e. number of 0.05s gaps in data
def shutter(band,t0,t1,tscale=1000.):
    return (str(baseURL)+
        'select shutter*0.05 from fGet'+str(band)+'Shutter('+
        str(long(t0*tscale))+','+str(long(t1*tscale))+')'+str(formatURL))

def shutdead(band,t0,t1,tscale=1000.):
	return (str(baseURL)+
        'SELECT shutter*0.05 FROM fGetNUVShutter('+str(long(t0*tscale))+','+
        str(long(t1*tscale))+
        ') AS time UNION ALL SELECT SUM(dt) * 0.0000057142857142857145 / ('+
        repr(t1)+'-'+repr(t0)+
        ') AS dead FROM(SELECT count(*) AS dt FROM NUVPhotonsNULLV'
        ' WHERE time BETWEEN '+str(long(t0*tscale))+' AND '+
        str(long(t1*tscale))+' UNION ALL SELECT count(*) AS dt'
        ' FROM NUVPhotonsV WHERE time BETWEEN '+str(long(t0*tscale))+' AND '+
        str(long(t1*tscale))+') x'+
        str(formatURL))

def exptime(band,t0,t1,stepsz=1.,tscale=1000.):
    return (str(baseURL)+
        'select * from fGet'+str(band)+'EffectiveExposureTime('+
        str(long(t0*tscale))+','+str(long(t1*tscale))+','+str(stepsz)+')'+
        str(formatURL))

def aspect(t0,t1,tscale=1000.):
    """Return the aspect information based on time range."""
    return (str(baseURL)+
        'select eclipse, filename, time, ra, dec, twist, flag, ra0, dec0,'
        ' twist0 from aspect where time between '+str(long(t0*tscale))+
        ' and '+str(long(t1*tscale))+' order by time'+str(formatURL))

# Return the aspect information based on eclipse
def aspect_ecl(eclipse):
    return (str(baseURL)+
        'select eclipse, filename, time, ra, dec, twist, flag, ra0, dec0,'
        ' twist0 from aspect where eclipse='+str(eclipse)+' order by time'+
        str(formatURL))

# Return data within a box centered on ra0, dec0 with sides of length 2*radius
def box(band,ra0,dec0,t0,t1,radius,tscale=1000.):
    return (str(baseURL)+
        'select time,ra,dec from '+str(band)+'PhotonsV where time between '+
        str(long(t0*tscale))+' and '+str(long(t1*tscale))+' and ra between '+
        repr(ra0-radius)+' and '+repr(ra0+radius)+' and dec between '+
        repr(dec0-radius)+' and '+repr(dec0+radius)+' and flag=0'+
        str(formatURL))

# Return data within a rectangle centered on ra0, dec0
def rect(band,ra0,dec0,t0,t1,ra,dec,tscale=1000.):
    return (str(baseURL)+
        'select time,ra,dec from fGetObjFromRect'+str(band)+'('+
        repr(ra0-ra/2.)+','+repr(ra0+ra/2.)+','+repr(dec0-dec/2.)+','+
        repr(dec0+dec/2.)+','+str(long(t0*tscale))+','+
        str(long(t1*tscale))+',0)'+str(formatURL))
