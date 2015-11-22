import numpy as np
import pandas as pd
import os
# gPhoton specific
import gQuery
from gQuery import tscale
import MCUtils as mc
import dbasetools as dbt # fGetTimeRanges(), compute_exptime()
import galextools as gxt # compute_flat_scale()
import cal

def gphot_params(band,skypos,radius,annulus=None,verbose=0.,detsize=1.25,
    stepsz=None,trange=None):
    """Populate a dict() with parameters that are constant over all bins."""
    return {'band':band,'ra0':skypos[0],'dec0':skypos[1],'skypos':skypos,
            'trange':trange,'radius':radius,'annulus':annulus,
            'stepsz':stepsz,'verbose':verbose,'detsize':detsize,
            'apcorrect1':gxt.apcorrect1(radius,band),
            'apcorrect2':gxt.apcorrect2(radius,band),
            'detbg':gxt.detbg(mc.area(radius),band)}

def xieta2colrow(xi, eta, band, detsize=1.25):
    """Convert detector xi, eta into col, row."""
    flat, flatinfo = cal.flat(band)
    # should be able to get npix from the header...
    npixx = flat.shape[0]
    npixy = flat.shape[1]
    pixsz = flatinfo['CDELT2']
    flatfill = detsize/(npixx*pixsz)
    col = ((( xi/36000.)/(detsize/2.)*flatfill + 1.)/2. * npixx)
    row = (((eta/36000.)/(detsize/2.)*flatfill + 1.)/2. * npixy)
    # You could theoretically drop a cut on detector position / detsize here...
    # Also, is this cut absolutely necessary? I think it's already been taken
    #  care of by the flag==0 assertion in the SQL query.
    #cut = ((col > 0.) & (col < flat.shape[0]-1) &
    #       (row > 0.) & (row < flat.shape[1]-1))
    #cut = np.where(ix == True)
    #ix = np.where((1.25/800.)*mc.distance(col,row,400,400)=detsize)
    return col, row

def hashresponse(band,events,verbose=0):
    """Given detector xi, eta, return the response at each position."""
    # Hash out the response correction
    if verbose:
        mc.print_inline("Applying the response correction.")
    flat, _ = cal.flat(band)
    events['col'], events['row'] = xieta2colrow(
                                            events['xi'], events['eta'], band)
    events['flat'] = flat[np.array(events['col'], dtype='int16'),
                          np.array(events['row'], dtype='int16')]
    events['scale'] = gxt.compute_flat_scale(events['t'], band)
    # TODO: Separately do the binlinearly interpolated response
    events['response'] = (events['flat']*events['scale'])
    return events

# Should this be moved to dbasetools?
def read_photons(photonfile,ra0,dec0,tranges,radius,verbose=0,
        colnames=['t','x','y','xa','ya','q','xi','eta','ra','dec','flags']):
    """Read a photon list file and return a python dict() with the expected
    format.
    """
    if verbose:
        print 'Reading photon list file: {f}'.format(f=photonfile)
    data = pd.io.parsers.read_csv(photonfile,names=colnames)
    ra,dec = np.array(data['ra']),np.array(data['dec'])
    angsep = mc.angularSeparation(ra0,dec0,ra,dec)
    ix = np.array([])
    for trange in tranges:
        if verbose:
            print trange
        cut = np.where((angsep<=radius) & (np.isfinite(angsep)))[0]
        ix = np.concatenate((ix,cut),axis=0)
    events = {'t':np.array(data['t'][ix])/tscale,
              'ra':np.array(data['ra'][ix]),
              'dec':np.array(data['dec'][ix]),
              'xi':np.array(data['xi'][ix]),
              'eta':np.array(data['eta'][ix]),
              'x':np.array(data['x'][ix]),
              'y':np.array(data['y'][ix])}
    return events

# This should be moved to dbasetools.
def query_photons(band,ra0,dec0,tranges,radius,verbose=0,flag=0):
    """Retrieve photons within an aperture from the database."""
    stream = []
    if verbose:
        print "Retrieving photons within {rad} degrees of [{r}, {d}]".format(
                                                        rad=radius,r=ra0,d=dec0)
    for trange in tranges:
        if verbose:
            mc.print_inline((" and between {t0} and {t1}.").format(
                                                    t0=trange[0],t1=trange[1]))
        thisstream = gQuery.getArray(
            gQuery.allphotons(band, ra0, dec0, trange[0], trange[1], radius,
                                        flag=flag), verbose=verbose,retries=100)
        stream.extend(thisstream)
    stream = np.array(stream, 'f8').T
    colnames = ['t', 'ra', 'dec', 'xi', 'eta', 'x', 'y']
    dtypes = ['f8', 'f8', 'f8', 'f4', 'f4', 'f4', 'f4']
    cols = map(np.asarray, stream, dtypes)
    events = dict(zip(colnames, cols))
    events['t']/=tscale # Adjust the timestamp by tscale
    return events

def pullphotons(band, ra0, dec0, tranges, radius, events={}, verbose=0,
                photonfile=None, flag=0):
    if photonfile:
        events = read_photons(photonfile, ra0, dec0, tranges, radius,
                              verbose=verbose)
    else:
        events = query_photons(band, ra0, dec0, tranges, radius,
                               verbose=verbose, flag=flag)

    events = hashresponse(band, events, verbose=verbose)
    return events

def aperture_error(counts,expt,bgcounts=0):
    return np.sqrt(counts/expt**2+bgcounts/expt**2)

def bg_sources(band,ra0,dec0,radius,margin=0.001):
    sources = gQuery.getArray(gQuery.mcat_sources(band,ra0,dec0,radius+margin,
                                                      maglimit=maskdepth))
    try:
        return {'ra':np.float32(np.array(sources)[:,0]),
                'dec':np.float32(np.array(sources)[:,1]),
                'fwhm':np.float32(np.array(sources)[:,7:9]),
                'radius':radius}
    except IndexError:
        return {'ra':np.array([]),'dec':np.array([]),
                'fwhm':np.array([]),'maglimit':maskdepth,'radius':radius}

def bg_mask_annulus(band,ra0,dec0,annulus,ras,decs,responses):
    ix = np.where((mc.angularSeparation(ra0,dec0,ras,decs)>=annulus[0]) &
                  (mc.angularSeparation(ra0,dec0,ras,decs)<=annulus[1]))
    return ras[ix],decs[ix],responses[ix]

def bg_mask_sources(band,ra0,dec0,ras,decs,responses,sources,maskradius=1.5):
    # At present, masks to 1.5 sigma where FWHM = 2.3548*sigma
    for i in range(len(sources['ra'])):
        ix = np.where(mc.angularSeparation(sources['ra'][i], sources['dec'][i],
                ras,decs)>=(maskradius/2.3548)*np.median(sources['fwhm'][i,:]))
        ras, decs, responses = ras[ix], decs[ix], responses[ix]
    return ras,decs,responses

def bg_mask(band,ra0,dec0,annulus,ras,decs,responses,sources):
    ras,decs,responses = bg_mask_annulus(band,ra0,dec0,annulus,ras,
                                         decs,responses)
    return bg_mask_sources(band,ra0,dec0,ras,decs,responses,sources)

def cheese_bg_area(band,ra0,dec0,annulus,sources,nsamples=10e5,ntests=10):
    # This is just a really naive Monte Carlo
    ratios = np.zeros(ntests)
    for i in range(ntests):
        ann_events = bg_mask_annulus(band,ra0,dec0,annulus,
             np.random.uniform(ra0-annulus[1],ra0+annulus[1],int(nsamples)),
             np.random.uniform(dec0-annulus[1],dec0+annulus[1],int(nsamples)),
             np.ones(nsamples))
        mask_events= bg_mask_sources(band,ra0,dec0,
                 ann_events[0],ann_events[1],ann_events[2],sources)
        try:
            ratios[i] = float(mask_events[2].sum())/float(ann_events[2].sum())
        except ZeroDivisionError:
            ratios[i] = 0.

    return (mc.area(annulus[1])-mc.area(annulus[0]))*ratios.mean()

# FIXME: This recomputes eff_area every pass at huge computational cost.
def cheese_bg(band,ra0,dec0,radius,annulus,ras,decs,responses,maskdepth=20.,
              maskradius=1.5,eff_area=False,sources=False):
    """ Returns the estimate number of counts (not count rate) within the
    aperture based upon a masked background annulus.
    """
    if not sources:
        sources = bg_sources(band,ra0,dec0,annulus[1],maskdepth=maskdepth)
    bg_counts = bg_mask(band,ra0,dec0,annulus,ras,decs,responses,
                        sources)[2].sum()
    if not eff_area:
        eff_area = cheese_bg_area(band,ra0,dec0,annulus,sources)
    return mc.area(radius)*bg_counts/eff_area if eff_area else 0.

def reduce_lcurve(bin_ix,region_ix,data,function,dtype='float64'):
    # Produces light curve columns by iteratively applying `function` to
    #  `data` within `region_ix` over `bin_ix`.
    bin_num = np.unique(bin_ix)
    output = np.empty(len(bin_num))
    for i,b in enumerate(bin_num):
        try:
            output[i] = function(data[np.where(bin_ix[region_ix]==b)])
        except ValueError:
            output[i] = np.nan
        except IndexError:
            output[i] = np.nan
        except:
            raise
    return np.array(output,dtype=dtype)

def maskwarning(band,bin_ix,events,verbose=0,mapkey='H'):
    # Test if any given events are _near_ a masked detector region.
    maps = {'H':cal.mask,'E':cal.flat}
    img,_ = maps[mapkey](band,buffer=True)
    ix = np.where(img[np.array(events['photons']['col'][bin_ix],dtype='int16'),
                np.array(events['photons']['row'][bin_ix],dtype='int16')]==0)
    return True if len(ix[0]) else False

def lowresponsewarning(bin_ix,events,verbose=0,ratio=0.7):
    ix = np.where(events['photons']['response'][bin_ix]<0.7)
    return True if len(ix[0]) else False

def exptimewarning(bin_ix,events,verbose=0,ratio=0.5):
    return (events['exptime'][bin_ix]/
                            (events['t1'][bin_ix]-events['t0'][bin_ix])<ratio)

def nonlinearitywarning(band,bin_ix,events,verbose=0):
    # Flags countrates above the 10% local nonlinearty dropoff per the
    # calibration paper.
    cps_10p_rolloff = {'NUV':311,'FUV':109}
    cps = events['flat_counts'][bin_ix]/events['exptime'][bin_ix]
    return True if cps>=cps_10p_rolloff[band] else False

def detedgewarning(bin_ix,events,verbose=0,valid_detrad=0.5):
    ix=np.where(mc.distance(events['photons']['col'][bin_ix],
        events['photons']['row'][bin_ix],400,400)*gxt.aper2deg(4)>=valid_detrad)
    return True if len(ix[0]) else False

def getflags(band,bin_ix,events,verbose=0):
    """Pass flags if data meets conditions that are likely to create
    misleading photometry. The flags are binary, with bins set as follows:
    1 - "hotspot" - events in pixels contiguous to a hotspot masked region
    2 - "mask edge" - events in pixels contiguous to the detector edge
    4 - "exptime" - bin contains < 50% exposure time coverage
    8 - "respose" - events weighted with response < 0.7
    16 - "nonlinearity" - local countrate exceeds 10% response dropoff
    32 - "detector edge" - events outside of 0.55 degrees of detector center

    """
    bin_num = np.unique(bin_ix)
    flags = np.zeros(len(bin_num))
    for i,b in enumerate(bin_num):
        try:
            ix = bin_ix[np.where(bin_ix==bin)]
            flags[i]+=1 if maskwarning(band,ix,events,mapkey='H',
                                                    verbose=verbose) else 0
            flags[i]+=2 if maskwarning(band,ix,events,mapkey='E',
                                                    verbose=verbose) else 0
            flags[i]+=4 if exptimewarning(i,events,verbose=verbose) else 0
            flags[i]+=8 if lowresponsewarning(ix,events,verbose=verbose) else 0
            flags[i]+=16 if nonlinearitywarning(band,i,events,
                                                    verbose=verbose) else 0
            flags[i]+=32 if detedgewarning(ix,events,verbose=verbose) else 0
        except IndexError:
            return np.array([np.nan])
        except:
            raise
    return np.array(flags)

def quickmag(band, ra0, dec0, tranges, radius, annulus=None, data={},
             stepsz=None, verbose=0,detsize=1.25,coadd=False, photonfile=None):
    if verbose:
        mc.print_inline("Retrieving all of the target events.")
    searchradius = radius if annulus is None else annulus[1]
    data = pullphotons(band, ra0, dec0, tranges, searchradius,
                       verbose=verbose)#, photonfile=photonfile)
    if verbose:
        mc.print_inline("Binning data according to requested depth.")
    # Multiple ways of defining bins
    try:
        trange = [np.array(tranges).min(),np.array(tranges).max()]
    except ValueError:
        trange = tranges
    if coadd:
        bins = np.array(trange)
    elif stepsz:
        bins = np.append(np.arange(trange[0], trange[1], stepsz), max(trange))
    else:
        bins = np.unique(np.array(tranges).flatten())

    lcurve = {'params':gphot_params(band,[ra0,dec0],radius,annulus=annulus,
                                    verbose=verbose,
                                    detsize=detsize,stepsz=stepsz,
                                    trange=trange)}

    """This is equivalent in function to np.digitize(data['t'],bins) except
    that it's much, much faster. See numpy issue #2656 at
    https://github.com/numpy/numpy/issues/2656
    """
    bin_ix = np.searchsorted(bins,data['t'],"right")
    try:
        lcurve['t0'] = bins[np.unique(bin_ix)-1]
        lcurve['t1'] = bins[np.unique(bin_ix)]
        lcurve['exptime'] = np.array(
            [dbt.compute_exptime(band,trange,skypos=[ra0,dec0],
                verbose=verbose,coadd=coadd)
                    for trange in zip(lcurve['t0'],lcurve['t1'])])
    except IndexError:
        if np.isnan(data['t']):
            print "No valid data available in {t}".format(t=tranges)
        lcurve['t0'] = np.array([np.nan])
        lcurve['t1'] = np.array([np.nan])
        lcurve['exptime'] = np.array([0])

    angSep = mc.angularSeparation(ra0, dec0, data['ra'], data['dec'])

    aper_ix = np.where(angSep <= radius)
    lcurve['t0_data']=reduce_lcurve(bin_ix,aper_ix,data['t'],np.min)
    lcurve['t1_data']=reduce_lcurve(bin_ix,aper_ix,data['t'],np.max)
    lcurve['t_mean']=reduce_lcurve(bin_ix,aper_ix,data['t'],np.mean)
    lcurve['counts']=reduce_lcurve(bin_ix,aper_ix,data['t'],len)
    lcurve['flat_counts']=reduce_lcurve(bin_ix,aper_ix,
                                                1./data['response'],np.sum)
    lcurve['responses']=reduce_lcurve(bin_ix,aper_ix,data['response'],np.mean)
    lcurve['detxs']=reduce_lcurve(bin_ix,aper_ix,data['col'],np.mean)
    lcurve['detys']=reduce_lcurve(bin_ix,aper_ix,data['row'],np.mean)
    lcurve['detrad'] = mc.distance(lcurve['detxs'],lcurve['detys'],400,400)
    lcurve['racent']=reduce_lcurve(bin_ix,aper_ix,data['ra'],np.mean)
    lcurve['deccent']=reduce_lcurve(bin_ix,aper_ix,data['dec'],np.mean)

    # FIXME: This just uses the average background over all visits until
    # the MCAT query gets fixed to return the time ranges.
    lcurve['mcat_bg']=lcurve['exptime']*(np.ones(len(lcurve['counts']))*
                        dbt.mcat_skybg(band,[ra0,dec0],radius,verbose=verbose))
    if annulus is not None:
        annu_ix = np.where((angSep > annulus[0]) & (angSep <= annulus[1]))
        lcurve['bg_counts'] = reduce_lcurve(bin_ix,annu_ix,data['t'],len)
        lcurve['bg_flat_counts']=reduce_lcurve(
                                    bin_ix,annu_ix,data['response'],np.sum)
        lcurve['bg'] = (mc.area(radius)*lcurve['bg_flat_counts'] /
                                    (mc.area(annulus[1])-mc.area(annulus[0])))
    else:
        lcurve['bg_counts'] = np.zeros(len(lcurve['counts']))
        lcurve['bg_flat_counts'] = np.zeros(len(lcurve['counts']))
        lcurve['bg'] = np.zeros(len(lcurve['counts']))

    lcurve['cps'] = lcurve['flat_counts']/lcurve['exptime']
    lcurve['cps_err'] = aperture_error(lcurve['flat_counts'],lcurve['exptime'])
    lcurve['cps_bgsub'] = (lcurve['flat_counts']-
                                            lcurve['bg'])/lcurve['exptime']
    lcurve['cps_bgsub_err'] = aperture_error(
        lcurve['flat_counts'],lcurve['exptime'],bgcounts=lcurve['bg'])
    lcurve['cps_mcatbgsub'] = (lcurve['flat_counts']-
                                        lcurve['mcat_bg'])/lcurve['exptime']
    lcurve['cps_mcatbgsub_err'] = aperture_error(
        lcurve['flat_counts'],lcurve['exptime'],bgcounts=lcurve['mcat_bg'])
    lcurve['flux'] = gxt.counts2flux(lcurve['cps'],band)
    lcurve['flux_err'] = gxt.counts2flux(lcurve['cps_err'],band)
    lcurve['flux_bgsub'] = gxt.counts2flux(lcurve['cps_bgsub'],band)
    lcurve['flux_bgsub_err'] = gxt.counts2flux(lcurve['cps_bgsub_err'],band)
    lcurve['flux_mcatbgsub'] = gxt.counts2flux(lcurve['cps_mcatbgsub'],band)
    lcurve['flux_mcatbgsub_err'] = gxt.counts2flux(
                                            lcurve['cps_mcatbgsub_err'],band)

    # NOTE: These conversions to mag can throw logarithm warnings if the
    # background is brighter than the source, resuling in a negative cps which
    # then gets propagated as a magnitude of NaN.
    lcurve['mag'] = gxt.counts2mag(lcurve['cps'],band)
    lcurve['mag_err_1'] = (lcurve['mag'] -
                        gxt.counts2mag(lcurve['cps'] + lcurve['cps_err'],band))
    lcurve['mag_err_2'] = (gxt.counts2mag(lcurve['cps'] -
                                    lcurve['cps_err'],band) - lcurve['mag'])
    lcurve['mag_bgsub'] = gxt.counts2mag(lcurve['cps_bgsub'],band)
    lcurve['mag_bgsub_err_1'] = (lcurve['mag_bgsub'] -
            gxt.counts2mag(lcurve['cps_bgsub'] + lcurve['cps_bgsub_err'],band))
    lcurve['mag_bgsub_err_2'] = (gxt.counts2mag(lcurve['cps_bgsub'] -
                            lcurve['cps_bgsub_err'],band) - lcurve['mag_bgsub'])
    lcurve['mag_mcatbgsub'] = gxt.counts2mag(lcurve['cps_mcatbgsub'],band)
    lcurve['mag_mcatbgsub_err_1'] = (lcurve['mag_mcatbgsub'] -
                                gxt.counts2mag(lcurve['cps_mcatbgsub'] +
                                            lcurve['cps_mcatbgsub_err'],band))
    lcurve['mag_mcatbgsub_err_2'] = (gxt.counts2mag(lcurve['cps_mcatbgsub'] -
                lcurve['cps_mcatbgsub_err'],band) - lcurve['mag_mcatbgsub'])

    lcurve['photons'] = data

    lcurve['flags'] = getflags(band,bin_ix,lcurve,verbose=verbose)

    return lcurve

def getcurve(band, ra0, dec0, radius, annulus=None, stepsz=None, lcurve={},
             trange=None, tranges=None, verbose=0, coadd=False, minexp=1.,
             maxgap=1.,photonfile=None, detsize=1.1):
    skyrange = [np.array(annulus).max().tolist() if annulus else radius,
                np.array(annulus).max().tolist() if annulus else radius,]
    if verbose:
        mc.print_inline("Getting exposure ranges.")
    if tranges is None:
        tranges = dbt.fGetTimeRanges(band, [ra0, dec0], trange=trange,
                maxgap=maxgap, minexp=minexp, verbose=verbose, detsize=detsize)
    elif not np.array(tranges).shape:
        print "No exposure time at this location: [{ra},{dec}]".format(
                                                            ra=ra0,dec=dec0)
    lcurve = quickmag(band, ra0, dec0, tranges, radius, annulus=annulus,
                              stepsz=stepsz, verbose=verbose, coadd=coadd)

    if verbose:
        mc.print_inline("Done.")
        mc.print_inline("")
    return lcurve

def write_curve(band, ra0, dec0, radius, csvfile=None, annulus=None,
                stepsz=None, trange=None, tranges=None, verbose=0, coadd=False,
                iocode='wb',detsize=1.1,overwrite=False,
                minexp=1.,maxgap=1.,photonfile=None):
    data = getcurve(band, ra0, dec0, radius, annulus=annulus, stepsz=stepsz,
                    trange=trange, tranges=tranges, verbose=verbose,
                    coadd=coadd, minexp=minexp, maxgap=maxgap,
                    photonfile=photonfile, detsize=detsize)
    if csvfile:
        columns = ['t0','t1','exptime','t_mean','t0_data','t1_data','cps',
                   'cps_err','cps_bgsub','counts','flat_counts','bg',
                   'mag','mag_err_1','mag_err_2','mag_bgsub','flux',
                   'flux_err','flux_bgsub','detx','dety','detrad','response',
                   'flags']
        try:
            test=pd.DataFrame({'t0':data['t0'],'t1':data['t1'],
                           't_mean':data['t_mean'],'t0_data':data['t0_data'],
                           't1_data':data['t1_data'],'exptime':data['exptime'],
                           'cps':data['cps'],'cps_err':data['cps_err'],
                           'cps_bgsub':data['cps_bgsub'],
                           'counts':data['counts'],
                           'flat_counts':data['flat_counts'],
                           'bg':data['bg'],'mag':data['mag'],
                           'mag_err_1':data['mag_err_1'],
                           'mag_err_2':data['mag_err_2'],
                           'mag_bgsub':data['mag_bgsub'],
                           'flux':data['flux'],
                           'flux_err':data['flux_err'],
                           'flux_bgsub':data['flux_bgsub'],
                           'detx':data['detxs'],'dety':data['detys'],
                           'detrad':data['detrad'],'response':data['responses'],
                           'flags':data['flags']
                           })
        except:
            raise
            if verbose>1:
                print 'Unable to build dataframe.'
        try:
            test.to_csv(csvfile,index=False,mode=iocode,columns=columns)
        except:
            print 'Did not write to: '+str(csvfile)
    else:
        if verbose>2:
            print "No CSV file requested."
        if verbose or (not verbose and not csvfile):
            print "AB Magnitudes:               "
            print data['mag']
    return data
