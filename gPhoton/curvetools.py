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

def gphot_params(band,skypos,radius,annulus=None,
                 verbose=0.,detsize=1.25,stepsz=None,
                 trange=None,maskdepth=None,maskradius=None):
    """Populate a dict() with parameters that are constant over all bins."""
    return {'band':band,'ra0':skypos[0],'dec0':skypos[1],'skypos':skypos,
            'trange':trange,'radius':radius,'annulus':annulus,
            'stepsz':stepsz,'verbose':verbose,
            'maskdepth':maskdepth,'maskradius':maskradius,'detsize':detsize,
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
            mc.print_inline(" and between "+str(trange[0])+" and "+
                            str(trange[1])+".")
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

def bg_sources(band,ra0,dec0,radius,maskdepth=20.0,maskradius=1.5,margin=0.001):
    sources = gQuery.getArray(gQuery.mcat_sources(band,ra0,dec0,radius+margin,
                                                      maglimit=maskdepth))
    try:
        return {'ra':np.float32(np.array(sources)[:,0]),
                'dec':np.float32(np.array(sources)[:,1]),
                'fwhm':np.float32(np.array(sources)[:,7:9]),
                'maskdepth':maskdepth,'maskradius':maskradius,
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
        except:
            raise
    return np.array(output,dtype=dtype)

def maskwarning(band,bin_ix,events,verbose=0,mapkey='H'):
    # Test if any given events are _near_ a masked detector region.
    maps = {'H':cal.mask,'E':cal.flat}
    img,_ = maps[mapkey](band,buffer=True)
    ix = np.where(img[np.array(events['col'][bin_ix],dtype='int16'),
                            np.array(events['row'][bin_ix],dtype='int16')]==0)
    return True if len(ix[0]) else False

def getflags(band,bin_ix,events,verbose=0):
    """Pass flags if data meets conditions that are likely to create
    misleading photometry.
    H = coincident with a masked hotspot.
    E = overlaps the edge of the detector.
    """
    bin_num = np.unique(bin_ix)
    output = np.empty(len(bin_num),dtype='string')
    for i,b in enumerate(bin_num):
        try:
            ix = bin_ix[np.where(bin_ix==bin)]
            output[i] = '{H}{E}'.format(
                H='H' if maskwarning(band,ix,events,mapkey='H',
                                                    verbose=verbose) else '',
                E='E' if maskwarning(band,ix,events,mapkey='E',
                                                    verbose=verbose) else '')
        except:
            raise
    return np.array(output)

def quickmag(band, ra0, dec0, tranges, radius, annulus=None, data={},
             stepsz=None, verbose=0, maskdepth=20.0,
             maskradius=1.5,detsize=1.25,coadd=False, photonfile=None):
    if verbose:
        mc.print_inline("Retrieving all of the target events.")
    searchradius = radius if annulus is None else annulus[1]
    data = pullphotons(band, ra0, dec0, tranges, searchradius,
                       verbose=verbose)#, photonfile=photonfile)
    if verbose:
        mc.print_inline("Binning data according to requested depth.")
    # Multiple ways of defining bins
    trange = [np.array(tranges).min(),np.array(tranges).max()]
    if coadd:
        bins = np.array(trange)
    elif stepsz:
        bins = np.append(np.arange(trange[0], trange[1], stepsz), max(trange))
    else:
        bins = np.unique(np.array(tranges).flatten())

    lcurve = {'params':gphot_params(band,[ra0,dec0],radius,annulus=annulus,
                                    verbose=verbose,
                                    detsize=detsize,stepsz=stepsz,
                                    trange=trange,maskdepth=maskdepth,
                                    maskradius=maskradius)}

    # This is equivalent in function to np.digitize(data['t'],bins) except
    # that it's much, much faster. See numpy issue #2656.
    bin_ix = np.searchsorted(bins,data['t'],"right")
    angSep = mc.angularSeparation(ra0, dec0, data['ra'], data['dec'])

    lcurve['t0'] = bins[np.unique(bin_ix)-1]
    lcurve['t1'] = bins[np.unique(bin_ix)]

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
    lcurve['racent']=reduce_lcurve(bin_ix,aper_ix,data['ra'],np.mean)
    lcurve['deccent']=reduce_lcurve(bin_ix,aper_ix,data['dec'],np.mean)

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

    lcurve['exptime'] = np.array(
        [dbt.compute_exptime(band,trange,skypos=[ra0,dec0],
            verbose=verbose,coadd=coadd)
                for trange in zip(lcurve['t0'],lcurve['t1'])])

    lcurve['flags'] = getflags(band,bin_ix,data,verbose=verbose)

    lcurve['photons'] = data
    return lcurve

def getcurve(band, ra0, dec0, radius, annulus=None, stepsz=None, lcurve={},
             trange=None, tranges=None, verbose=0, coadd=False, minexp=1.,
             maxgap=1., maskdepth=20, maskradius=1.5,
             photonfile=None, detsize=1.1):
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
    # FIXME: Everything goes to hell if no exposure time is available...
    # TODO: Add an ability to specify or exclude specific time ranges
    if verbose:
        mc.print_inline("Moving to photon level operations.")
    # FIXME: This error handling is hideous.
    try:
        lcurve = quickmag(band, ra0, dec0, tranges, radius, annulus=annulus,
                          stepsz=stepsz, verbose=verbose, coadd=coadd)
        lcurve['cps'] = lcurve['flat_counts']/lcurve['exptime']
        lcurve['cps_err'] = np.sqrt(lcurve['flat_counts'])/lcurve['exptime']
        lcurve['cps_bgsub'] = (lcurve['flat_counts']-
                               lcurve['bg'])/lcurve['exptime']
        lcurve['mag'] = gxt.counts2mag(lcurve['cps'],band)
        lcurve['mag_err_1'] = (lcurve['mag'] - gxt.counts2mag(lcurve['cps'] + lcurve['cps_err'],band))
        lcurve['mag_err_2'] = gxt.counts2mag(lcurve['cps'] - lcurve['cps_err'],band) - lcurve['mag']
        lcurve['mag_bgsub'] = gxt.counts2mag(lcurve['cps_bgsub'],band)
        lcurve['flux'] = gxt.counts2flux(lcurve['cps'],band)
        lcurve['flux_err'] = gxt.counts2flux(lcurve['cps_err'],band)
        lcurve['flux_bgsub'] = gxt.counts2flux(lcurve['cps_bgsub'],band)
        lcurve['detrad'] = mc.distance(lcurve['detxs'],lcurve['detys'],400,400)
    except ValueError:
        lcurve['cps']=[]
        lcurve['cps_err']=[]
        lcurve['cps_bgsub']=[]
        lcurve['mag']=[]
        lcurve['mag_err_1']=[]
        lcurve['mag_err_2']=[]
        lcurve['mag_bgsub']=[]
        lcurve['flux']=[]
        lcurve['flux_err']=[]
        lcurve['flux_bgsub']=[]
        lcurve['detrad']=[]
    if verbose:
        mc.print_inline("Done.")
        mc.print_inline("")
    return lcurve

def write_curve(band, ra0, dec0, radius, csvfile=None, annulus=None,
                stepsz=None, trange=None, tranges=None, verbose=0, coadd=False,
                iocode='wb',detsize=1.1,overwrite=False,
                minexp=1.,maxgap=1.,maskdepth=20.,maskradius=1.5,
                photonfile=None):
    data = getcurve(band, ra0, dec0, radius, annulus=annulus, stepsz=stepsz,
                    trange=trange, tranges=tranges, verbose=verbose,
                    coadd=coadd, minexp=minexp, maxgap=maxgap,
                    maskdepth=maskdepth, maskradius=maskradius,
                    photonfile=photonfile, detsize=detsize)
    if csvfile:
        columns = ['t0','t1','exptime','t_mean','t0_data','t1_data','cps',
                   'cps_err','cps_bgsub','counts','flat_counts','bg',
                   'mag','mag_err_1','mag_err_2','mag_bgsub','flux',
                   'flux_err','flux_bgusb','detx','dety','detrad','response',
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
