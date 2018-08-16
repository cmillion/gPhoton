"""
.. module:: gQuery
   :synopsis: Defines and constructs common queries that are passed to the
       GALEX databases (esp: photon, aspect, and MCAT) at MAST.
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
from builtins import str
# gPhoton imports.
import gPhoton.CalUtils as CalUtils
from gPhoton.MCUtils import manage_requests2
from gPhoton.galextools import isPostCSP
from gPhoton import time_id

# ------------------------------------------------------------------------------
# To save space, times in the database are "integer-ized" by multiplying by 1000
tscale = 1000.

# The following three global variables are used in constructing a properly
# formatted query to the MAST database. Don't change them unless you know what
# you're doing!
baseURL = ('https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/Galex'
           'PhotonListQueryTest?query=')
baseDB = 'GPFCore.dbo'
MCATDB = 'GR6Plus7.dbo'

# All queries from the same _run_ of the photon tools should have identical
# time_id, providing a quick way to troubleshoot issues on the server side.
formatURL = ' -- '+str(time_id)+'&format=extjs'
# ------------------------------------------------------------------------------

# The photon even timestamps are stored in the database at the precision level
#  of SQL's BIGINT. This truncated (not rounded) some timestamps at the level
#  of 1ms. Most timestamps have a resolution of only 5ms except for rare high
#  resolution visits, and even in that case the extra precision does not
#  matter for science. To make gAperture consistent with the database, we'll
#  truncate times at 1ms for queries.
def truncate(n):
    return str(n*tscale).split('.')[0]

# ------------------------------------------------------------------------------
def hasNaN(query):
    """
    Check if there is NaN in a query (or any string) and, if so, raise an
        exception because that probably indicates that something has gone wrong.

    :param query: The query string to check.

    :type query: str
    """

    if 'NaN' in query:
        raise RuntimeError("Malformed query: contains NaN values.")

    return
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def getValue(query, verbose=0, retries=100):
    """
    Manage a database call which returns a single value.

    :param query: The query to run.

    :type query: str

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int

    :returns: requests.Response or None -- The response from the server. If the
        query does not receive a response, returns None.
    """

    hasNaN(query)

    out = manage_requests2(query, maxcnt=retries, verbose=verbose)

    if out is not None:
        try:
            out = float(out.json()['data']['Tables'][0]['Rows'][0][0])
        except ValueError:
            try:
                out = str(out.json()['data']['Tables'][0]['Rows'][0][0])
            except:
                print('Failed: {q}'.format(q=query))
                raise
        except:
            print('Failed: {q}'.format(q=query))
            raise
        return out
    else:
        print('Failed: {q}'.format(q=query))
        raise ValueError("Query never finished on server, run with verbose"
                         " turned on for more info.")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def getArray(query, verbose=0, retries=100):
    """
    Manage a database call which returns an array of values.

    :param query: The query to run.

    :type query: str

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int

    :param retries: Number of query retries to attempt before giving up.

    :type retries: int

    :returns: requests.Response or None -- The response from the server. If the
        query does not receive a response, returns None.
    """

    hasNaN(query)

    out = manage_requests2(query, maxcnt=retries, verbose=verbose)

    if out is not None:
        try:
            out = out.json()['data']['Tables'][0]['Rows']
        except:
            print('Failed: {q}'.format(q=query))
            raise
        return out
    else:
        print('Failed: {q}'.format(q=query))
        raise ValueError("Query never finished on server, run with verbose"
                         " turned on for more info.")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def mcat_sources(band, ra0, dec0, radius, maglimit=20.):
    """
    Return the MCAT coadd sources given sky position and search radius
        (and optional lower magnitude limit).
        Columns are:
        [0,RA],[1,Dec],[2,NUV_mag],[3,FUV_mag],[4,FoV_radius],[5,NUV_skybg],
        [6,FUV_skybg],[7,NUV_FWHM_world],[8,FUV_FWHM_world],
        [9:15,FUV_mag_aper_1:7],[16:22,NUV_mag_aper_1:7]
        [23:29,FUV_magerr_aper_1:7],[30:36,NUV_magerr_aper1:7]

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: The right ascension, in degrees, around which to search.

    :type ra0: float

    :param dec0: The declination, in degrees, around which to search.

    :type dec0: float

    :param radius: The radius within which to search for MCAT sources, in
        degrees.

    :type radius: float

    :param maglimit: The NUV faint limit to return MCAT sources for.

    :type maglimit: float

    :returns: str -- The query to submit to the database.
    """

    # 1=nuv, 2=fuv, 3=both
    bandflag = 1 if band == 'NUV' else 2

    # fGetNearbyObjEq takes radius in arcminutes
    # [Future]: Add exposure time.
    return (
        str(baseURL)+
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
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def obstype(objid):
    """
    Get the dither pattern type based on the object id.

    :param objid: The MCAT Object ID to return the observation type data from.

    :type objid: long

    :returns: str -- The query to submit to the database.
    """

    return (
        '{baseURL}select distinct vpe.mpstype as survey, vpe.tilename,'
        ' vpe.photoextractid, vpe.petal, vpe.nlegs, vpe.leg, vpe.eclipse,'
        ' vpe.img, vpe.subvis from {MCATDB}.visitPhotoextract as vpe inner'
        ' join {MCATDB}.imgrun as iv on vpe.photoextractid=iv.imgrunid'
        ' inner join {MCATDB}.visitphotoobjall as p on vpe.photoextractid'
        '=p.photoextractid where p.objid={objid}{formatURL}'.format(
            baseURL=baseURL, MCATDB=MCATDB, objid=objid, formatURL=formatURL))
# ------------------------------------------------------------------------------

def obstype_from_t(t):
    """
    Get the dither pattern type based on the time stamp.
    """
    return ("{baseURL}SELECT * from {baseDB}.fGetLegObsType({t})"
            "{formatURL}").format(baseURL=baseURL, baseDB=baseDB,
                                  t=truncate(t), formatURL=formatURL)

# -_----------------------------------------------------------------------------
def mcat_visit_sources(ra0, dec0, radius):
    """
    Return the MCAT per-visit sources given sky position and search radius.
        The columns are as follows:
        [0,objid],[1,ra],[2,dec],[3,NUV_mag],[4,FUV_mag],[5,FoV_radius],
        [6,NUV_skybg],[7,FUV_skybg],[8,NUV_FWHM],[9,FUV_FWHM],[10,FUV_expt],
        [11,NUV_expt],[12:18,FUV_mag_aper_1:7],[19:25,NUV_mag_aper_1:7],
        [26:32,FUV_magerr_aper_1:7],[33:39,NUV_magerr_aper_1:7],[40,Nobssecs],
        [41,Fobssecs],[42,NUV_artifact],[43,FUV_artifact],[44,FUV_obstart],
        [45,FUV_obsend],[46,NUV_obstart],[47,NUV_obsend],
        [48,FUV_ALPHA_J2000],[49,FUV_DELTA_J2000],
        [50,NUV_ALPHA_J2000],[51,NUV_DELTA_J2000]

    :param ra0: The right ascension, in degrees, around which to search.

    :type ra0: float

    :param dec0: The declination, in degrees, around which to search.

    :type dec0: float

    :param radius: The radius within which to search for MCAT sources, in
        degrees.

    :type radius: float

    :returns: str -- The query to submit to the database.
    """

    return (
        "{baseURL}select vpo.objid, ra, dec, nuv_mag, fuv_mag, fov_radius,"
        " nuv_skybg, fuv_skybg, nuv_fwhm_world, fuv_fwhm_world,"
        " vpe.fexptime, vpe.nexptime, fuv_mag_aper_1, fuv_mag_aper_2,"
        " fuv_mag_aper_3, fuv_mag_aper_4, fuv_mag_aper_5, fuv_mag_aper_6,"
        " fuv_mag_aper_7, nuv_mag_aper_1, nuv_mag_aper_2, nuv_mag_aper_3,"
        " nuv_mag_aper_4, nuv_mag_aper_5, nuv_mag_aper_6, nuv_mag_aper_7,"
        " fuv_magerr_aper_1, fuv_magerr_aper_2, fuv_magerr_aper_3,"
        " fuv_magerr_aper_4, fuv_magerr_aper_5, fuv_magerr_aper_6,"
        " fuv_magerr_aper_7, nuv_magerr_aper_1, nuv_magerr_aper_2,"
        " nuv_magerr_aper_3, nuv_magerr_aper_4, nuv_magerr_aper_5,"
        " nuv_magerr_aper_6, nuv_magerr_aper_7, nobssecs, fobssecs,"
        " nuv_artifact, fuv_artifact, vpe.fexpstar, vpe.fexpend,"
        " vpe.nexpstar, vpe.nexpend, FUV_ALPHA_J2000, FUV_DELTA_J2000,"
        " NUV_ALPHA_J2000, NUV_DELTA_J2000"
        " from {MCATDB}.visitphotoobjall as vpo"
        " inner join {MCATDB}.visitphotoextract as vpe on"
        " vpo.photoextractid=vpe.photoextractid inner join"
        " {MCATDB}.fGetNearbyVisitObjEq({ra0},{dec0},{radius}) as nb on"
        " vpo.objid=nb.objid inner join {MCATDB}.imgrun as i on"
        " vpe.photoExtractID=i.imgRunID{formatURL}".format(
            baseURL=baseURL, MCATDB=MCATDB, ra0=float(ra0), dec0=float(dec0),
            radius=radius*60., formatURL=formatURL))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def mcat_objid_search(objid):
    """
    Return a bunch of observation data for a visit level objid (ggoid).
        Doing the same for coadd level data is not yet supported.

    :param objid: The MCAT Object ID to return the observation data from.

    :type objid: long

    :returns: str -- The query to submit to the database.
    """

    return (
        str(baseURL)+'select objid, minPhotoObsDate, maxPhotoObsDate, obs_date,'
        ' obsdatim, nobssecs, fobssecs, nexptime, fexptime, nexpstar, nexpend,'
        ' fexpstar, fexpend from '+str(MCATDB)+'.visitphotoobjall as vp inner'
        ' join '+str(MCATDB)+'.imgrun as ir on vp.photoextractid=ir.imgrunid'
        ' inner join '+str(MCATDB)+'.visitphotoextract as vpe on'
        ' vp.photoextractid=vpe.photoextractid where objid = '+
        str(int(objid))+str(formatURL))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def exposure_ranges(band, ra0, dec0, t0=1, t1=10000000000000, detsize=1.25,
                    epsilon=0.001):
    """
    Returns a list of times (in one second increments) where data exists
        with an aspect solution within detsize of [ra0,dec0].

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: The right ascension, in degrees, around which to search.

    :type ra0: float

    :param dec0: The declination, in degrees, around which to search.

    :type dec0: float

    :param t0: The minimum time stamp to search for exposure ranges.

    :type t0: long

    :param t1: The maximum time stamp to search for exposure ranges.

    :type t1: long

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :param epsilon: Buffer on t1 to avoid missing the end value in the search.

    :type epsilon: float

    :returns: str -- The query to submit to the database.
    """

    # If band is set to False (or equivalent), search on both bands
    if not band:
        band = 'FUV/NUV'

    return (
        str(baseURL)+
        'select distinct time from '+str(baseDB)+
        '.fGetNearbyAspectEq('+repr(float(ra0))+','+repr(float(dec0))+',(('+
        str(detsize)+'/2.0)*60.0),'+
        truncate(t0)+','+truncate(t1+epsilon)+')'
        ' where band=\''+str(band)+'\' or band=\'FUV/NUV\' order by time'+
        str(formatURL))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def exposure_range(band, ra0, dec0, t0=1, t1=10000000000000):
    """
    Find time ranges for which data exists at a given position.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: The right ascension, in degrees, around which to search.

    :type ra0: float

    :param dec0: The declination, in degrees, around which to search.

    :type dec0: float

    :param t0: The minimum time stamp to search for exposure ranges.

    :type t0: long

    :param t1: The maximum time stamp to search for exposure ranges.

    :type t1: long

    :returns: str -- The query to submit to the database.
    """

    return (str(baseURL)+'select startTimeRange, endTimeRange'
            ' from '+str(baseDB)+
            '.fGetTimeRanges('+str(int(t0))+','+str(int(t1))+','+
            repr(float(ra0))+','+repr(float(dec0))+') where band=\''+
            str(band)+'\''+str(formatURL))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def aperture(band, ra0, dec0, t0, t1, radius):
    """
    Integrate counts over an aperture at a position.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: The right ascension, in degrees, around which to search.

    :type ra0: float

    :param dec0: The declination, in degrees, around which to search.

    :type dec0: float

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :param radius: The radius within which to integrate counts, in
        degrees.

    :type radius: float

    :returns: str -- The query to submit to the database.
    """

    return (str(baseURL)+
            'Select  sum(photonCount) from '+str(baseDB)+
            '.fGetNearbyObjEqCount'+str(band)+'('+repr(float(ra0))+','+
            repr(float(dec0))+','+str(radius)+','+truncate(t0)+','+
            truncate(t1)+',0)'+str(formatURL))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def deadtime1(band, t0, t1, flag=False):
    """
    Return the global counts of non-NULL data.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :param flag: If true, return only flag=0 data, else return all non-NULL.

    :type flag: bool

    :returns: str -- The query to submit to the database.
    """

    return (
        '{baseURL}select count(*) from {baseDB}.{band}PhotonsV where '+
        'time >= {t0} and time < {t1}'+
        '{flag}{formatURL}').format(
            baseURL=baseURL, baseDB=baseDB, band=band, t0=truncate(t0),
            t1=truncate(t1), flag=' and flag=0' if flag else '',
            formatURL=formatURL)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def deadtime2(band, t0, t1):
    """
    Return the global counts for NULL data.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :returns: str -- The query to submit to the database.
    """

    return ('{baseURL}select count(*) from {baseDB}.{band}PhotonsNULLV where '+
            'time >= {t0} and time < {t1}{formatURL}').format(
                baseURL=baseURL, baseDB=baseDB, band=band,
                t0=truncate(t0), t1=truncate(t1),
                formatURL=formatURL)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def deadtime(band, t0, t1, feeclkratio=0.966, tec2fdead=5.52e-6):
    """
    Return the emperically determined deadtime correction based upon the
        global count rate.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :param feeclkratio: Ratio of Front End Electronics clock rates.

    :type feeclkratio: float

    :param tec2fdead: The nominal amount of time following an event that the
        detector is unable to detect another event.

    :type tec2fdead: float

    :returns: str -- The query to submit to the database.
    """

    scale = tec2fdead/feeclkratio

    return (str(baseURL)+
            'select sum(dt)*'+'%0.30f'%scale+' / ('+repr(t1)+'-'+repr(t0)+')'
            ' from(select count(*) as dt from '+str(baseDB)+'.'+str(band)+
            'PhotonsNULLV where time >= '+truncate(t0)+' and time < '+
            truncate(t1)+' union all select count(*) as dt from '+
            str(baseDB)+'.'+str(band)+
            'PhotonsV where time >= '+truncate(t0)+' and time < '+
            truncate(t1)+') x'+str(formatURL))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def globalcounts(band, t0, t1, flag=False):
    """
    Return the total number of detector events within a time range.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :param flag: If true, return only flag=0 data. Else return all non-NULL.

    :type flag: bool

    :returns: str -- The query to submit to the database.
    """

    return ('{baseURL}select count(t) from '+
            '(select time as t from {baseDB}.{band}PhotonsV '+
            'where time >= {t0} and time < {t1}{flag} union all '+
            'select time as t from {baseDB}.{band}PhotonsNULLV where time >= '+
            '{t0} and time < {t1}) x{formatURL}').format(
                baseURL=baseURL, baseDB=baseDB, band=band,
                t0=truncate(t0), t1=truncate(t1),
                flag=' and flag=0' if flag else '', formatURL=formatURL)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def alltimes(band, t0, t1):
    """
    Return the time stamps of every detector event within a time range.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :returns: str -- The query to submit to the database.
    """

    return ('{baseURL}select t from (select time as t from'
            ' {baseDB}.{band}PhotonsV '+
            'where time >= {t0} and time < {t1} union all select time as t '+
            'from {baseDB}.{band}PhotonsNULLV where time >= '+
            '{t0} and time < {t1}) x{formatURL}').format(
                baseURL=baseURL, baseDB=baseDB, band=band,
                t0=truncate(t0), t1=truncate(t1),
                formatURL=formatURL)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def uniquetimes(band, t0, t1, flag=False, null=False):
    """
    Return the unique timestamps for events within trange.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :param flag: If true, only return flag=0 data. Else return all non-NULL.

    :type flag: bool

    :param null: If true, query the null table.

    :returns: str -- The query to submit to the database.
    """

    if not null:
        return ('{baseURL}select time from {baseDB}.{band}ShutterPOTimeV '+
                'where time >= {t0} and time < {t1} order by '+
                'time{formatURL}').format(
                    baseURL=baseURL, baseDB=baseDB, band=band,
                    t0=truncate(t0), t1=truncate(t1),
                    formatURL=formatURL)
    else:
        return ('{baseURL}select distinct time from '+
                '{baseDB}.{band}PhotonsNULLV where '+
                'time >= {t0} and time < {t1} order by time'+
                '{formatURL}').format(
                    baseURL=baseURL, baseDB=baseDB, band=band,
                    t0=truncate(t0), t1=truncate(t1),
                    formatURL=formatURL)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def boxcount(band, t0, t1, xr, yr):
    """
    Find the number of events inside of a box defined by [xy] range in
        detector space coordinates. This is useful for pulling out stim events.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :param xr: The minimum and maximum x-values that define the box.

    :type xr: list

    :param yr: The minimum and maximum y-values that define the box.

    :type yr: list

    :returns: str -- The query to submit to the database.
    """

    return (str(baseURL)+'select count(*) from '+str(baseDB)+'.'+str(band)+
            'PhotonsNULLV where time >= '+truncate(t0)+' and time < '+
            truncate(t1)+' and x >= '+str(xr[0])+' and x < '+str(xr[1])+
            ' and y >= '+str(yr[0])+' and y < '+str(yr[1])+str(formatURL))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def detbox(band, t0, t1, xr, yr):
    """
    Return all events inside a box defined in detector space by [xy]
        range. Created as a sanity check for stim events.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :param xr: The minimum and maximum x-values that define the box.

    :type xr: list

    :param yr: The minimum and maximum y-values that define the box.

    :type yr: list

    :returns: str -- The query to submit to the database.
    """

    return ('{baseURL}select time,x,y,ya from {baseDB}.{band}PhotonsNULLV '+
            'where time >= {t0} and time < {t1} and '+
            'x >= {xmin} and x < {xmax} and y >= {ymin} and y < {ymax}'+
            '{formatURL}').format(
                baseURL=baseURL, baseDB=baseDB, band=band,
                t0=truncate(t0), t1=truncate(t1), xmin=xr[0],
                xmax=xr[1], ymin=yr[0], ymax=yr[1], formatURL=formatURL)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def stimcount(band, t0, t1, margin=[90.01, 90.01], aspum=68.754932/1000.,
              eclipse=None, null=True):
    """
    Return stim counts.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :param margin: X and Y lengths, in arcseconds, of a box within which
        to search for stim events.

    :type margin: list

    :param aspum: Arcseconds per micrometer (of detector).

    :type aspum: float

    :param eclipse: The eclipse to return stim counts for.

    :type eclipse: int

    :param null: If true, query the NULL table.

    :type null: bool

    :returns: str -- The query to submit to the database.
    """

    # [Future]: Convert t0 to eclipse

    if not eclipse:
        eclipse = 55000 if isPostCSP(t0) else 30000

    if isPostCSP(t0):
        margin[1] = 180.02

    avgstim = CalUtils.avg_stimpos(band, eclipse)

    return ('{baseURL}select count(*) from {baseDB}.{band}Photons{N}V '+
            'where time >= {t0} and time < {t1} and ('+
            '((x >= {x10} and x < {x11}) and (y >= {y10} and y < {y11})) or '+
            '((x >= {x20} and x < {x21}) and (y >= {y20} and y < {y21})) or '+
            '((x >= {x30} and x < {x31}) and (y >= {y30} and y < {y31})) or '+
            '((x >= {x40} and x < {x41}) and (y >= {y40} and y < {y41}))'+
            '){formatURL}').format(baseURL=baseURL, baseDB=baseDB, band=band,
                                   N='NULL' if null else '',
                                   t0=truncate(t0),
                                   t1=truncate(t1),
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
                                   y41=(avgstim['y4']+margin[1])/aspum,
                                   formatURL=formatURL)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def stimtimes(band, t0, t1, margin=[90.01, 90.01], aspum=68.754932/1000.,
              eclipse=None):
    """
    Return stim counts.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :param margin: X and Y lengths, in arcseconds, of a box in which to
        search for stim values.

    :type margin: list

    :param aspum: Arcseconds per micrometer (on detector).

    :type aspum: float

    :param eclipse: The eclipse to return stim counts for.

    :type eclipse: int

    :returns: str -- The query to submit to the database.
    """

    if not eclipse:
        eclipse = 55000 if isPostCSP(t0) else 30000

    if isPostCSP(t0):
        margin[1] = 180.02

    avgstim = CalUtils.avg_stimpos(band, eclipse)

    return ('{baseURL}select time from {baseDB}.{band}PhotonsNULLV '+
            'where time >= {t0} and time < {t1} and ('+
            '((x >= {x10} and x < {x11}) and (y >= {y10} and y < {y11})) or '+
            '((x >= {x20} and x < {x21}) and (y >= {y20} and y < {y21})) or '+
            '((x >= {x30} and x < {x31}) and (y >= {y30} and y < {y31})) or '+
            '((x >= {x40} and x < {x41}) and (y >= {y40} and y < {y41}))'+
            '){formatURL}').format(baseURL=baseURL, baseDB=baseDB, band=band,
                                   t0=truncate(t0),
                                   t1=truncate(t1),
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
                                   y41=(avgstim['y4']+margin[1])/aspum,
                                   formatURL=formatURL)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def boxcentroid(band, t0, t1, xr, yr):
    """
    Find the mean position of events inside of a box in detector space.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :param xr: The minimum and maximum x-values that define the box.

    :type xr: list

    :param yr: The minimum and maximum y-values that define the box.

    :type yr: list

    :returns: str -- The query to submit to the database.
    """

    return (str(baseURL)+'select avg(x), avg(y) from '+str(baseDB)+
            '.'+str(band)+'PhotonsNULLV where time >= '+
            truncate(t0)+' and time < '+truncate(t1)+
            ' and x >= '+str(xr[0])+' and x < '+str(xr[1])+' and y >= '+
            str(yr[0])+' and y < '+str(yr[1])+str(formatURL))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def boxtimes(band, t0, t1, xr, yr):
    """
    Get the list of times for events inside of a box in detector space.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :param xr: The minimum and maximum x-values that define the box.

    :type xr: list

    :param yr: The minimum and maximum y-values that define the box.

    :type yr: list

    :returns: str -- The query to submit to the database.
    """

    return (str(baseURL)+
            'select time from '+str(baseDB)+'.'+str(band)+
            'PhotonsNULLV where time >= '+
            truncate(t0)+' and time < '+truncate(t1)+
            ' and x >= '+
            str(xr[0])+' and x < '+str(xr[1])+' and y >= '+str(yr[0])+
            ' and y < '+str(yr[1])+str(formatURL))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def allphotons(band, ra0, dec0, t0, t1, radius, flag=0):
    """
    Grab the major columns for all events within an aperture.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: The right ascension, in degrees, around which to search.

    :type ra0: float

    :param dec0: The declination, in degrees, around which to search.

    :type dec0: float

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :param radius: The radius within which to search, in degrees.

    :type radius: float

    :param flag: Only return times with this flag value. Zero is nominal.
        NOTE: 'Flag' is not a reliable way to parse data at this time. You
        should compare event timestamps against the aspect file.

    :type flag: int

    :returns: str -- The query to submit to the database.
    """

    return ('{baseURL}select time,ra,dec,xi,eta,x,y,flag,q from '
            '{baseDB}.fGetNearbyObjEq{band}AllColumns({ra0},{dec0},{radius},'
            '{t0},{t1},{flag}){formatURL}').format(
                baseURL=baseURL, baseDB=baseDB, band=band, ra0=repr(float(ra0)),
                dec0=repr(float(dec0)), radius=radius, t0=truncate(t0),
                t1=truncate(t1), flag=flag, formatURL=formatURL)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def shutter(band, t0, t1):
    """
    Get shutter correction, i.e., number of 0.05-sec gaps in data.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :returns: str -- The query to submit to the database.
    """

    return (str(baseURL)+'select shutter*0.05 from '+str(baseDB)+'.fGet'+
            str(band)+'Shutter('+truncate(t0)+','+truncate(t1)+
            ')'+str(formatURL))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def aspect(t0, t1):
    """
    Return the aspect information based on a time range.

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :returns: str -- The query to submit to the database.
    """

    return (str(baseURL)+
            'select eclipse, filename, time, ra, dec, twist, flag, ra0, dec0,'
            ' twist0 from aspect where time >= '+truncate(t0)+
            ' and time < '+truncate(t1)+' order by time'+str(formatURL))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def aspect_ecl(eclipse):
    """
    Return the aspect information based upon an eclipse number.

    :param eclipse: The eclipse to return aspect information for.

    :type eclipse: int

    :returns: str -- The query to submit to the database.
    """

    return (str(baseURL)+
            'select eclipse, filename, time, ra, dec, twist, flag, ra0, dec0,'
            ' twist0 from aspect where eclipse='+str(eclipse)+' order by time'+
            str(formatURL))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def aspect_skypos(ra, dec, detsize=1.25):
    """
    Return the aspect information based upon sky position and det radius.

    :param ra: The right ascension to search on, in degrees.

    :type ra: float

    :param dec: The declination to search on, in degrees.

    :type dec: float

    :param detsize: Effective diameter, in degrees, of the field-of-view.

    :type detsize: float

    :returns: str -- The query to submit to the database.
    """

    return (str(baseURL)+
            "select eclipse, filename, time, ra, dec, twist, flag, ra0, dec0,"
            " twist0 from aspect where ra >= "+repr(ra-detsize/2.)+
            " and ra < "+repr(ra+detsize/2.)+" and dec >= "+
            repr(dec-detsize/2.)+" and dec < "+repr(dec+detsize/2.)+
            " order by time"+str(formatURL))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def box(band, ra0, dec0, t0, t1, radius, flag=0):
    """
    Return data within a box centered on ra0, dec0 with sides of length
        2*radius.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: The right ascension, in degrees, around which to search.

    :type ra0: float

    :param dec0: The declination, in degrees, around which to search.

    :type dec0: float

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :param radius: The radius within which to search, in degrees.

    :type radius: float

    :param flag: If true, only return flag=0 data. Else return all non-NULL.

    :type flag: int

    :returns: str -- The query to submit to the database.
    """

    # [Future]: Deprecate this and rename it skybox().

    return (str(baseURL)+'select time,ra,dec from '+str(baseDB)+'.'+str(band)+
            'PhotonsV where time >= '+
            truncate(t0)+' and time < '+truncate(t1)+
            ' and ra >= '+
            repr(ra0-radius)+' and ra < '+repr(ra0+radius)+' and dec >= '+
            repr(dec0-radius)+' and dec < '+repr(dec0+radius)+' and flag='+
            str(flag)+str(formatURL))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def skyrect(band, ra0, dec0, t0, t1, ra, dec, flag=0):
    """
    Extract photon data falling within a specified rectangle on the sky.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param ra0: The right ascension, in degrees, around which to search.

    :type ra0: float

    :param dec0: The declination, in degrees, around which to search.

    :type dec0: float

    :param t0: The minimum time stamp to search.

    :type t0: long

    :param t1: The maximum time stamp to search.

    :type t1: long

    :param ra: The length in degrees along RA describing the region of
        interest.

    :type ra: float

    :param dec: The length in degrees along Dec describing the region of
        interest.

    :type dec: float

    :param flag: Flag value of non-NULL data to return. Zero is nominal.

    :type flag: int

    :returns: str -- The query to submit to the database.
    """

    return ('{baseURL}select time,ra,dec,xi,eta,x,y from'
            ' {baseDB}.fGetObjFromRect{band}AllColumns({ra_min},{ra_max},'
            '{dec_min},{dec_max},{t0},{t1},{flag}){formatURL}'.format(
                baseURL=baseURL, baseDB=baseDB, band=band,
                ra_min=repr(ra0-ra/2.), ra_max=repr(ra0+ra/2.),
                dec_min=repr(dec0-dec/2.), dec_max=repr(dec0+dec/2.),
                t0=truncate(t0), t1=truncate(t1), flag=flag,
                formatURL=formatURL))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def raw_data_paths(eclipse):
    """
    Construct a query that returns a data structure containing the download
    paths

    :param eclipse: GALEX eclipse number.

    :type flag: int

    :returns: str -- The query to submit to the database.
    """
    return 'https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=spGetRawUrls {ecl}&format=extjs'.format(ecl=int(eclipse))
# ------------------------------------------------------------------------------
