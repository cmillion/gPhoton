import numpy as np
import pandas as pd
import os
# gPhoton specific
import gQuery
import MCUtils as mc
import dbasetools as dbt # fGetTimeRanges(), compute_exptime()
import galextools as gxt # compute_flat_scale()
import imagetools as it # sigmaclip_bg()
from FileUtils import flat_filename

def gphot_params(band,skypos,radius,annulus=None,calpath='../cal/',
                 verbose=0.,detsize=1.25,stepsz=None,
                 trange=None,maskdepth=None,maskradius=None):
    """Populate a dict() with parameters that are constant over all bins."""
    return {'band':band,'ra0':skypos[0],'dec0':skypos[1],'skypos':skypos,
            'trange':trange,'radius':radius,'annulus':annulus,
            'stepsz':stepsz,'calpath':calpath,'verbose':verbose,
            'maskdepth':maskdepth,'maskradius':maskradius,'detsize':detsize,
            'apcorrect1':gxt.apcorrect1(radius,band),
            'apcorrect2':gxt.apcorrect2(radius,band),
            'detbg':gxt.detbg(mc.area(radius),band)}

def xieta2colrow(xi, eta, calfile, detsize=1.25):
    """Convert detector xi, eta into col, row."""
    flat = mc.get_fits_data(calfile)
    flatinfo = mc.get_fits_header(calfile)
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
    return col, row

def hashresponse(band,events,calpath='../cal/',verbose=0):
    """Given detector xi, eta, return the response at each position."""
    # Hash out the response correction
    if verbose:
        mc.print_inline("Applying the response correction.")
    flat = mc.get_fits_data(flat_filename(band, calpath))
    events['col'], events['row'] = xieta2colrow(events['xi'], events['eta'],
                                                flat_filename(band, calpath))
    events['flat'] = flat[np.array(events['col'], dtype='int16'),
                          np.array(events['row'], dtype='int16')]
    events['scale'] = gxt.compute_flat_scale(events['t'], band)
    # TODO: Separately do the binlinearly interpolated response
    events['response'] = (events['flat']*events['scale'])
    return events

def pullphotons(band, ra0, dec0, tranges, radius, events={}, verbose=0,
                tscale=1000., calpath='../cal/', chunksz=10e6):
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
            gQuery.allphotons(band, ra0, dec0, trange[0], trange[1], radius),
                              verbose=verbose,retries=100)
#        if not stream:
#            continue
        stream.extend(thisstream)
    stream = np.array(stream, 'f8').T
    colnames = ['t', 'ra', 'dec', 'xi', 'eta']
    dtypes = ['f8', 'f8', 'f8', 'f4', 'f4']
    cols = map(np.asarray, stream, dtypes)
    events = dict(zip(colnames, cols))
    events['t']/=tscale # Adjust the timestamp by tscale
    events = hashresponse(band, events, calpath=calpath, verbose=verbose)
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
    #mc.print_inline('Estimating area of masked background annulus.')
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
    #mc.print_inline('Swiss cheesing the background annulus.')
    if not sources:
        sources = bg_sources(band,ra0,dec0,annulus[1],maskdepth=maskdepth)
    bg_counts = bg_mask(band,ra0,dec0,annulus,ras,decs,responses,
                        sources)[2].sum()
    #mc.print_inline('Numerically integrating area of masked annulus.')
    if not eff_area:
        eff_area = cheese_bg_area(band,ra0,dec0,annulus,sources)
    return mc.area(radius)*bg_counts/eff_area if eff_area else 0.

def quickmag(band, ra0, dec0, tranges, radius, annulus=None, data={},
             stepsz=None, calpath='../cal/', verbose=0, maskdepth=20.0,
             maskradius=1.5,detsize=1.25,coadd=False,sigmaclip=3.):
    if verbose:
        mc.print_inline("Retrieving all of the target events.")
    trange = [np.array(tranges).min(),np.array(tranges).max()]
    try:
        searchradius = annulus[1]
    except TypeError:
        searchradius = radius
    data = pullphotons(band, ra0, dec0, tranges, searchradius,
                       verbose=verbose, calpath=calpath)
    if verbose:
        mc.print_inline("Isolating source from background.")
    angSep = mc.angularSeparation(ra0, dec0, data['ra'], data['dec'])
    if verbose:
        mc.print_inline("Binning data according to requested depth.")
    # Multiple ways of defining bins
    if coadd:
        bins = np.array(trange)
    elif stepsz:
        bins = np.append(np.arange(min(trange), max(trange), stepsz),
                                                                max(trange))
    else:
        bins = np.unique(np.array(tranges).flatten())
    # This is equivalent in function to np.digitize(data['t'],bins) except
    # that it's much, much faster. See numpy issue #2656.
    ix = np.searchsorted(bins,data['t'],"right")
    # Initialize histogrammed arrays
    # FIXME: allocate these from a dict of constructors
    lcurve_cols = ['counts', 'sources', 'bg_counts','responses',
                   'detxs', 'detys', 't0_data', 't1_data', 't_mean', 'racent',
                   'deccent']
    lcurve = {'params':gphot_params(band,[ra0,dec0],radius,annulus=annulus,
                                    calpath=calpath,verbose=verbose,
                                    detsize=detsize,stepsz=stepsz,
                                    trange=trange,maskdepth=maskdepth,
                                    maskradius=maskradius)}
    for col in lcurve_cols:
        lcurve[col] = np.zeros(len(bins)-1)
    # FIXME: Bottleneck. There's probably a way to do this without looping.
    # Don't bother looping through anything with no data.
    lcurve['bg'] = {'simple':np.zeros(len(bins)-1),
                    'cheese':np.zeros(len(bins)-1),
                    'sigmaclip':np.zeros(len(bins)-1)}
    if annulus is not None:
        lcurve['bg']['sources'] = bg_sources(band,ra0,dec0,annulus[1],
                                             maskdepth=maskdepth)
        lcurve['bg']['eff_area'] = cheese_bg_area(band,ra0,dec0,annulus,
                                                  lcurve['bg']['sources'])
    else:
        lcurve['bg']['sources'] = None
        lcurve['bg']['eff_area'] = 0.
    if verbose:
        mc.print_inline("Populating histograms.")
    for cnt,i in enumerate(np.unique(ix)):
        # Exclude data outside of the bins in searchsorted.
        if i-1<0 or i==len(bins):
            continue
        if verbose:
            mc.print_inline('Binning {i} of {l}.'.format(
                                                    i=cnt,l=len(np.unique(ix))))
        t_ix = np.where(ix==i)
        # TODO: Optionally limit data to specific parts of detector.
        rad_ix = np.where((angSep <= radius) & (ix == i))
        lcurve['t0_data'][i-1] = data['t'][t_ix and rad_ix].min()
        lcurve['t1_data'][i-1] = data['t'][t_ix and rad_ix].max()
        lcurve['t_mean'][i-1] = data['t'][t_ix and rad_ix].mean()
        lcurve['counts'][i-1] = len(rad_ix[0])
        lcurve['sources'][i-1] = (1./data['response'][rad_ix]).sum()
        lcurve['responses'][i-1] = data['response'][rad_ix].mean()
        lcurve['detxs'][i-1] = data['col'][rad_ix].mean()
        lcurve['detys'][i-1] = data['row'][rad_ix].mean()
        lcurve['racent'][i-1] = data['ra'][rad_ix].mean()
        lcurve['deccent'][i-1] = data['dec'][rad_ix].mean()
        if annulus is not None:
            ann_ix = np.where((angSep > annulus[0]) &
                              (angSep <= annulus[1]) & (ix == i))
            lcurve['bg_counts'][i-1] = len(ann_ix[0])
            # Background is reported as counts within the aperture
            lcurve['bg']['simple'][i-1] = (mc.area(radius) *
                (1./data['response'][ann_ix]).sum() /
                (mc.area(annulus[1])-mc.area(annulus[0])))
            lcurve['bg']['cheese'][i-1] = cheese_bg(band, ra0, dec0, radius,
                annulus, data['ra'][t_ix], data['dec'][t_ix],
                data['response'][t_ix], maskdepth=maskdepth,
                eff_area=lcurve['bg']['eff_area'],
                sources=lcurve['bg']['sources'])
            lcurve['bg']['sigmaclip'][i-1] = it.sigmaclip_bg(data,radius,
                                            annulus,[ra0,dec0],sigmaclip=1)
        else:
            lcurve['bg_counts'][i-1]=0.
            lcurve['bg']['simple'][i-1]=0.
            lcurve['bg']['cheese'][i-1]=0.
            lcurve['bg']['sigmaclip'][i-1]=0.
    # Only return bins that contain data.
    ix = np.where((np.isfinite(lcurve['sources'])) &
                  (np.array(lcurve['sources']) > 0))
    lcurve['t0'] = bins[ix]
    lcurve['t1'] = bins[ix[0]+1]
    for col in lcurve_cols:
        lcurve[col] = lcurve[col][ix]
    if annulus is not None:
        lcurve['bg']['simple']=lcurve['bg']['simple'][ix]
        lcurve['bg']['cheese']=lcurve['bg']['cheese'][ix]
        lcurve['bg']['sigmaclip']=lcurve['bg']['sigmaclip'][ix]
    else:
        lcurve['bg']['simple']=0.
        lcurve['bg']['cheese']=0.
        lcurve['bg']['sigmaclip']=0.

    lcurve['exptime'] = np.array(
        [dbt.compute_exptime(band,trange,skypos=[ra0,dec0],
                             verbose=verbose,coadd=coadd)
            for trange in zip(lcurve['t0_data'],lcurve['t1_data'])])
    if verbose:
        mc.print_inline("Returning curve data.")
    lcurve['photons'] = data
    return lcurve

def getcurve(band, ra0, dec0, radius, annulus=None, stepsz=None, lcurve={},
             trange=None, tranges=None, verbose=0, coadd=False, minexp=1.,
             maxgap=1., calpath='../cal/', maskdepth=20, maskradius=1.5,
             sigmaclip=3.):
    if verbose:
        mc.print_inline("Getting exposure ranges.")
    if tranges is None:
        tranges = dbt.fGetTimeRanges(band, [ra0, dec0], trange=trange,
                                 maxgap=maxgap, minexp=minexp, verbose=verbose)
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
                          stepsz=stepsz, verbose=verbose, coadd=coadd,
                          calpath=calpath, maskdepth=maskdepth,
                          maskradius=maskradius,sigmaclip=sigmaclip)
        lcurve['cps'] = lcurve['sources']/lcurve['exptime']
        lcurve['cps_bgsub'] = (lcurve['sources']-
                               lcurve['bg']['simple'])/lcurve['exptime']
        lcurve['cps_bgsub_cheese'] = (lcurve['sources']-
                               lcurve['bg']['cheese'])/lcurve['exptime']
        lcurve['cps_bgsub_sigmaclip'] = (lcurve['sources']-
                            lcurve['bg']['sigmaclip'])/lcurve['exptime']
        lcurve['mag'] = gxt.counts2mag(lcurve['cps'],band)
        lcurve['mag_bgsub'] = gxt.counts2mag(lcurve['cps_bgsub'],band)
        lcurve['mag_bgsub_cheese'] = gxt.counts2mag(
                                            lcurve['cps_bgsub_cheese'],band)
        lcurve['mag_bgsub_sigmaclip'] = gxt.counts2mag(
                                        lcurve['cps_bgsub_sigmaclip'],band)
        lcurve['flux'] = gxt.counts2flux(lcurve['cps'],band)
        lcurve['flux_bgsub'] = gxt.counts2flux(lcurve['cps_bgsub'],band)
        lcurve['flux_bgsub_cheese'] = gxt.counts2flux(
                                            lcurve['cps_bgsub_cheese'],band)
        lcurve['flux_bgsub_sigmaclip'] = gxt.counts2flux(
                                        lcurve['cps_bgsub_sigmaclip'],band)
        lcurve['detrad'] = mc.distance(lcurve['detxs'],lcurve['detys'],400,400)

    except ValueError:
        lcurve['cps']=[]
        lcurve['cps_bgsub']=[]
        lcurve['cps_bgsub_cheese']=[]
        lcurve['cps_bgsub_sigmaclip']=[]
        lcurve['mag']=[]
        lcurve['mag_bgsub']=[]
        lcurve['mag_bgsub_cheese']=[]
        lcurve['mag_bgsub_sigmaclip']=[]
        lcurve['flux']=[]
        lcurve['flux_bgsub']=[]
        lcurve['flux_bgsub_cheese']=[]
        lcurve['flux_bgusb_sigmaclip']=[]
        lcurve['detrad']=[]
    if verbose:
        mc.print_inline("Done.")
        mc.print_inline("")
    return lcurve

def write_curve(band, ra0, dec0, radius, csvfile=None, annulus=None,
                stepsz=None, trange=None, tranges=None, verbose=0, coadd=False,
                iocode='wb',calpath='../cal/',detsize=1.25,clobber=False,
                minexp=1.,maxgap=1.,maskdepth=20.,maskradius=1.5,sigmaclip=3.):
    data = getcurve(band, ra0, dec0, radius, annulus=annulus, stepsz=stepsz,
                    trange=trange, tranges=tranges, verbose=verbose,
                    coadd=coadd, minexp=minexp, maxgap=maxgap, calpath=calpath,
                    maskdepth=maskdepth, maskradius=maskradius,
                    sigmaclip=sigmaclip)
    print tranges
    if csvfile:
        columns = ['t0','t1','exptime','mag_bgsub_cheese','t_mean','t0_data',
                   't1_data','cps','counts','bg','mag','mag_bgsub',
                   'flux','flux_bgsub','flux_bgsub_cheese','bg_cheese']
        test=pd.DataFrame({'t0':data['t0'],'t1':data['t1'],
                           't_mean':data['t_mean'],'t0_data':data['t0_data'],
                           't1_data':data['t1_data'],'exptime':data['exptime'],
                           'cps':data['cps'],'counts':data['counts'],
                           'bg':data['bg']['simple'],'mag':data['mag'],
                           'mag_bgsub':data['mag_bgsub'],
                           'mag_bgsub_cheese':data['mag_bgsub_cheese'],
                           'flux':data['flux'],
                           'flux_bgsub':data['flux_bgsub'],
                           'flux_bgsub_cheese':data['flux_bgsub_cheese'],
                           'bg_cheese':data['bg']['cheese']})#,
                          #'bg_sigmaclip':data['bg']['sigmaclip'],
                          #'cps_bgsub_sigmaclip':data['cps_bgsub_sigmaclip'],
                          #'mag_bgsub_sigmaclip':data['mag_bgsub_sigmaclip'],
                          #'flux_bgsub_sigmaclip':data['flux_bgsub_sigmaclip']})
        try:
            test.to_csv(csvfile,index=False,mode=iocode,columns=columns)
        except:
            print 'Failed to write to: '+str(csvfile)
    else:
        if verbose>2:
            print "No CSV file requested."
        if verbose or (not verbose and not csvfile):
            print "AB Magnitudes:               "
            print data['mag']
    return data
