import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import gQuery
import MCUtils as mc
import dbasetools as dbt
import galextools as gxt
from FileUtils import flat_filename

def gphot_params(band,skypos,radius,annulus=[False,False],calpath='../cal/',verbose=False,detsize=1.25,stepsz=False,trange=[1,1000000000000]):
    return {'band':band,'skypos':skypos,'radius':radius,'annulus':annulus,'calpath':calpath,'verbose':verbose,'detsize':detsize,'stepsz':stepsz,'trange':trange}

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
    # You could theoretically drop a cut on detector position in here...
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

def pullphotons(band, ra0, dec0, t0, t1, radius, events={}, verbose=0,
                tscale=1000.,calpath='../cal/'):
    """Retrieve photons within an aperture from the database."""
    if verbose:
        mc.print_inline("Retrieving all photons at ["+str(ra0)+", "+
                        str(dec0)+"] within a radius of "+str(radius)+
                        " and between "+str(t0)+" and "+str(t1)+".")
    stream = gQuery.getArray(
        gQuery.allphotons(band, ra0, dec0, t0, t1, radius), verbose=verbose)
    if not stream:
        return events
    events['t'] = np.array(stream,dtype='float64')[:,0]/tscale
    # The float64 precision _is_ significant for RA / Dec.
    events['ra'] = np.array(stream,dtype='float64')[:,1]
    events['dec'] = np.array(stream,dtype='float64')[:,2]
    events['xi'] = np.array(stream,dtype='float32')[:,3]
    events['eta']= np.array(stream,dtype='float32')[:,4]
    events = hashresponse(band, events, calpath=calpath, verbose=verbose)
    return events

def bg_sources(band,ra0,dec0,radius,maglimit=28.):
    sources = gQuery.getArray(gQuery.mcat_sources(band,ra0,dec0,radius,
                                                      maglimit=maglimit))
    return {'ra':np.float32(np.array(sources)[:,0]),
            'dec':np.float32(np.array(sources)[:,1]),
            'fwhm':np.float32(np.array(sources)[:,7:9]),
            'maglimit':maglimit,'radius':radius}

def bg_mask_annulus(band,ra0,dec0,annulus,ras,decs,responses):
    ix = np.where((mc.angularSeparation(ra0,dec0,ras,decs)>=annulus[0]) &
                  (mc.angularSeparation(ra0,dec0,ras,decs)<=annulus[1]))
    return ras[ix],decs[ix],responses[ix]

def bg_mask_sources(band,ra0,dec0,ras,decs,responses,sources):
    for i in range(len(sources['ra'])):
        ix = np.where(mc.angularSeparation(sources['ra'][i], sources['dec'][i],
                      ras,decs)>=sources['fwhm'][i,:].max())
        ras, decs, responses = ras[ix], decs[ix], responses[ix]
    return ras,decs,responses

def bg_mask(band,ra0,dec0,annulus,ras,decs,responses,sources):
    ras,decs,responses = bg_mask_annulus(band,ra0,dec0,annulus,ras,
                                         decs,responses)
    return bg_mask_sources(band,ra0,dec0,ras,decs,responses,sources)

def cheese_bg_area(band,ra0,dec0,annulus,sources,nsamples=10e4,ntests=10):
    mc.print_inline('Estimating area of masked background annulus.')
    ratios = np.zeros(ntests)
    #print band,ra0,dec0,annulus,sources.keys()
    for i in range(ntests):
        #print i
        ann_events = bg_mask_annulus(band,ra0,dec0,annulus,
                 np.random.uniform(ra0-annulus[1],ra0+annulus[1],nsamples),
                 np.random.uniform(dec0-annulus[1],dec0+annulus[1],nsamples),
                 np.ones(nsamples))
        mask_events= bg_mask_sources(band,ra0,dec0,
                 ann_events[0],ann_events[1],ann_events[2],sources)
        ratios[i] = float(mask_events[2].sum())/float(ann_events[2].sum())
    return (mc.area(annulus[1])-mc.area(annulus[0]))*ratios.mean()

def cheese_bg(band,ra0,dec0,radius,annulus,ras,decs,responses,maglimit=28.,
              eff_area=False,sources=False):
    """ Returns the estimate number of counts (not count rate) within the
    aperture based upon a masked background annulus.
    """
    mc.print_inline('Swiss cheesing the background annulus.')
    if not sources:
        sources = bg_sources(band,ra0,dec0,annulus[1],maglimit=maglimit)
    bg_counts = bg_mask(band,ra0,dec0,annulus,ras,decs,responses,
                        sources)[2].sum()
    mc.print_inline('Numerically integrating area of masked annulus.')
    if not eff_area:
        eff_area = cheese_bg_area(band,ra0,dec0,annulus,sources)
    return mc.area(radius)*bg_counts/eff_area

def quickmag(band, ra0, dec0, trange, radius, annulus=False, data={},
             stepsz=False, calpath='../cal/', verbose=0, maglimit=28.):
    if verbose:
        mc.print_inline("Retrieving all of the target events.")
    data = pullphotons(band, ra0, dec0, trange[0], trange[1],
                       radius if not annulus else annulus[1])
    if verbose:
        mc.print_inline("Isolating source from background.")
    angSep = mc.angularSeparation(ra0, dec0, data['ra'], data['dec'])
    if verbose:
        mc.print_inline("Binning data according to requested depth.")
    bins = (np.append(np.arange(min(trange), max(trange), stepsz), max(trange)) 
            if stepsz else np.array([min(trange), max(trange)]))
    ix = np.digitize(data['t'],bins)
    # Initialize histogramed arrays.
    # FIXME: allocate these from a dict of constructors
    lcurve_cols = ['counts', 'sources', 'bg_counts', 'bkgrnds', 'responses',
                   'detxs', 'detys', 't0_data', 't1_data', 't_mean']
    lcurve = {'params':{'band':band,'ra0':ra0,'dec0':dec0,'trange':trange,
              'radius':radius,'annulus':annulus,'stepsz':stepsz,
              'calpath':calpath,'verbose':verbose,'maglimit':maglimit}}
    for col in lcurve_cols:
        lcurve[col] = np.zeros(len(bins)-1)
    # FIXME: Bottleneck. There's probably a way to do this without looping.
    # Don't bother looping through anything with no data.
    lcurve['bg'] = {'simple':np.zeros(len(bins)-1),
                    'cheese':np.zeros(len(bins)-1)}
    lcurve['bg']['sources'] = bg_sources(band,ra0,dec0,annulus[1],maglimit=maglimit)
    lcurve['bg']['eff_area'] = cheese_bg_area(band,ra0,dec0,annulus,lcurve['bg']['sources'])
    for cnt,i in enumerate(np.unique(ix)):
        if verbose:
            mc.print_inline('Binning '+str(i)+' of '+str(len(ix))+'.')
        t_ix = np.where(ix==i)
        lcurve['t0_data'][i-1] = data['t'][t_ix].min()
        lcurve['t1_data'][i-1] = data['t'][t_ix].max()
        lcurve['t_mean'][i-1] = data['t'][t_ix].mean()
        # TODO: Optionally limit data to specific parts of detector.
        rad_ix = np.where((angSep <= radius) & (ix == i))
        lcurve['counts'][i-1] = len(rad_ix[0])
        lcurve['sources'][i-1] = (1./data['response'][rad_ix]).sum()
        lcurve['responses'][i-1] = data['response'][rad_ix].mean()
        lcurve['detxs'][i-1] = data['col'][rad_ix].mean()
        lcurve['detys'][i-1] = data['row'][rad_ix].mean()
        ann_ix = np.where((angSep > annulus[0]) &
                          (angSep <= annulus[1]) & (ix == i))
        lcurve['bg_counts'][i-1] = len(ann_ix[0])
        lcurve['bg']['simple'][i-1] = mc.area(radius)*(1./data['response'][ann_ix]).sum()/(mc.area(annulus[1])-mc.area(annulus[0]))
        lcurve['bg']['cheese'][i-1] = cheese_bg(band, ra0, dec0, radius,
            annulus, data['ra'][t_ix], data['dec'][t_ix],
            data['response'][t_ix], maglimit=maglimit,
            eff_area=lcurve['bg']['eff_area'],sources=lcurve['bg']['sources'])
    # Only return bins that contain data.
    ix = np.where((np.isfinite(lcurve['sources'])) &
                  (np.array(lcurve['sources']) > 0))
    lcurve['t0'] = bins[ix]
    lcurve['t1'] = bins[ix[0]+1]
    for col in lcurve_cols:
        lcurve[col] = lcurve[col][ix]
    lcurve['bg']['simple']=lcurve['bg']['simple'][ix]
    lcurve['bg']['cheese']=lcurve['bg']['cheese'][ix]
    if verbose:
        mc.print_inline("Returning curve data.")
    # Return: t0, t1, counts, background counts, events,
    # background events, mean response, mean detector x, mean detector y
    lcurve['photons'] = data
    return lcurve

def getcurve(band, ra0, dec0, radius, annulus=False, stepsz=False, lcurve={},
             trange=[1,1000000000000], verbose=0):
    if verbose:
        mc.print_inline("Getting exposure ranges.")
    tranges = dbt.fGetTimeRanges(band, [ra0, dec0], trange=trange, maxgap=100,
                                 minexp=100, verbose=verbose)
    # TODO: Add an ability to specify or exclude specific time ranges
    if verbose:
        mc.print_inline("Moving to photon level operations.")
    lcurve = quickmag(band, ra0, dec0, [tranges.min(), tranges.max()], radius,
                      annulus=annulus, stepsz=stepsz, verbose=verbose)
    if verbose:
        mc.print_inline("Scaling background to aperture.")
    # FIXME: This exposure time calculation is a bottleneck!
    # The overhead is in the actual network call, so you could construct
    # a single large query that returns an array with all of the exposure
    # time parameters (for all bins) at once.
    if verbose:
        mc.print_inline("Computing exposure times... could take a while...")
    lcurve['exptime'] = np.array(
        [dbt.compute_exptime(band, [lcurve['t0_data'][i], lcurve['t1_data'][i]])#, skypos=[ra0,dec0])
         for i in range(len(lcurve['t0_data']))])
    #print lcurve['counts'].shape,lcurve['bg']['cheese'].shape,lcurve['exptime'].shape
    lcurve['cps'] = (lcurve['counts']-lcurve['bg']['cheese'])/lcurve['exptime']
    if verbose:
        mc.print_inline("Done.")
        mc.print_inline("")
    return lcurve

def new_write_curve(band, ra0, dec0, radius, outfile, annulus=False,
                    stepsz=False, trange=[1,1000000000000], verbose=0,
                    iocode='wb',calpath='../cal/',detsize=1.25):
    data = getcurve(band, ra0, dec0, radius, annulus=annulus, stepsz=stepsz,
                    trange=trange, verbose=verbose)
    cols = ['counts', 'sources', 'bg_counts', 'bkgrnds', 'responses',
            'detxs', 'detys', 't0_data', 't1_data', 't_mean']
    test=pd.DataFrame({'t0':data['t0'],'t1':data['t1'],'t_mean':data['t_mean'],'t0_data':data['t0_data'],'t1_data':data['t1_data'],'exptime':data['exptime'],'cps':data['cps'],'counts':data['counts'],'bg':data['bg']['cheese']})
    # TODO: clobber check
    try:
        test.to_csv(outfile,index=False)
    except:
        print 'Failed to write to: '+str(outfile)
    return
