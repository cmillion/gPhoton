import numpy as np
import matplotlib.pyplot as plt

import gQuery
import MCUtils as mc
import dbasetools as dbt
import galextools as gxt
from FileUtils import flat_filename

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

def pullphotons(band, ra0, dec0, t0, t1, radius, data={}, verbose=0,
                tscale=1000.):
    """Retrieve photons within an aperture from the database."""
    if verbose:
        mc.print_inline("Retrieving all photons at ["+str(ra0)+", "+
                        str(dec0)+"] within a radius of "+str(radius)+
                        " and between "+str(t0)+" and "+str(t1)+".")
    stream = gQuery.getArray(
        gQuery.allphotons(band, ra0, dec0, t0, t1, radius), verbose=verbose)
    if not stream:
        return data
    data['t'] = np.array(stream,dtype='float64')[:,0]/tscale
    # The float64 precision _is_ significant for RA / Dec.
    data['ra'] = np.array(stream,dtype='float64')[:,1]
    data['dec'] = np.array(stream,dtype='float64')[:,2]
    data['xi'] = np.array(stream,dtype='float32')[:,3]
    data['eta']= np.array(stream,dtype='float32')[:,4]
    return data

def hashresponse(band,xi,eta,data={},calpath='../cal/',verbose=0):
    """Given detector xi, eta, return the response at each position."""
    # Hash out the response correction
    if verbose:
        mc.print_inline("Applying the response correction.")
    flat = mc.get_fits_data(flat_filename(band, calpath))
    # TODO: Do you really need ix here?
    data['col'], data['row'] = xieta2colrow(data['xi'], data['eta'],
    flat_filename(band, calpath))
    data['flat'] = flat[np.array(data['col'], dtype='int16'),
    np.array(data['row'], dtype='int16')]
    data['scale'] = gxt.compute_flat_scale(data['t'], band)
    data['response'] = (data['flat']*data['scale'])
    return data

def quickmag(band, ra0, dec0, trange, radius, annulus=False, data={},
             stepsz=False, calpath='../cal/', verbose=0):
    if verbose:
        mc.print_inline("Retrieving all of the target events.")
    data = pullphotons(band, ra0, dec0, trange[0], trange[1],
                       radius if not annulus else annulus[1], data=data)
    # Hash out the response correction
    data = hashresponse(band, data['xi'], data['eta'], data=data,
                        calpath=calpath, verbose=verbose)
    if verbose:
        mc.print_inline("Isolating source from background.")
    angSep = mc.angularSeparation(ra0, dec0, data['ra'], data['dec'])
    if verbose:
        mc.print_inline("Binning data according to requested depth.")
    bins = (np.append(np.arange(min(trange), max(trange), stepsz), max(trange)) 
            if stepsz else np.array([min(trange), max(trange)]))
    ix = np.digitize(data['t'],bins)
    # Initialize histogramed arrays.
    counts = np.zeros(len(bins)-1) # Number of events in aperture
    sources = np.zeros(len(bins)-1) # Scaled integral of events in aperture
    bg_counts = np.zeros(len(bins)-1) # Number of events in annulus
    bkgrnds = np.zeros(len(bins)-1) # Scaled integral of events in annulus
    responses = np.zeros(len(bins)-1) # Mean response in aperture
    detxs = np.zeros(len(bins)-1) # Mean detector x position in aperture
    detys = np.zeros(len(bins)-1) # Mean detector y position in aperture
    t0_data = np.zeros(len(bins)-1) # Minimum time value of data
    t1_data = np.zeros(len(bins)-1) # Maximum time value of data
    # FIXME: Bottleneck. There's probably a way to do this without looping.
    # Don't bother looping through anything with no data.
    for cnt,i in enumerate(np.unique(ix)):
        if verbose:
            mc.print_inline('Binning '+str(i)+' of '+str(len(ix))+'.')
        t_ix = np.where(ix==i)
        t0_data[i-1] = data['t'][t_ix].min()
        t1_data[i-1] = data['t'][t_ix].max()
        # TODO: Optionally limit data to specific parts of detector.
        rad_ix = np.where((angSep <= radius) & (ix == i))
        counts[i-1] = len(rad_ix[0])
        sources[i-1] = (1./data['response'][rad_ix]).sum()
        responses[i-1] = data['response'][rad_ix].mean()
        detxs[i-1] = data['col'][rad_ix].mean()
        detys[i-1] = data['row'][rad_ix].mean()
        # TODO: Exclude regions around nearby sources.
        ann_ix = np.where((angSep > annulus[0]) &
                          (angSep <= annulus[1]) & (ix == i))
        bg_counts[i-1] = len(ann_ix[0])
        bkgrnds[i-1] = (1./data['response'][ann_ix]).sum()
    # Only return bins that contain data.
    ix = np.where((np.isfinite(sources)) & (np.array(sources) > 0))
    if verbose:
        mc.print_inline("Returning curve data.")
    # Return: t0, t1, counts, background counts, events,
    # background events, mean response, mean detector x, mean detector y
    return (bins[ix], bins[ix[0]+1], t0_data[ix], t1_data[ix], sources[ix],
            bkgrnds[ix], counts[ix], bg_counts[ix], responses[ix], detxs[ix],
            detys[ix], data)

def getcurve(band, ra0, dec0, radius, annulus=False, stepsz=False, curve={},
             trange=[1,1000000000000], verbose=0):
    if verbose:
        mc.print_inline("Getting exposure ranges.")
    tranges = dbt.fGetTimeRanges(band, [ra0, dec0], trange=trange, maxgap=100,
                                 minexp=100, verbose=verbose)
    # TODO: Add an ability to specify or exclude specific time ranges
    if verbose:
        mc.print_inline("Moving to photon level operations.")
    (curve['t0'], curve['t1'], curve['t0_data'], curve['t1_data'], 
            curve['counts'], curve['bg'], curve['events'],
            curve['bg_events'], curve['response'], curve['detx'],
            curve['dety'], curve['data']) = quickmag(band, ra0, dec0,
            [tranges.min(), tranges.max()], radius,
            annulus=annulus, stepsz=stepsz, verbose=verbose)
    if verbose:
        mc.print_inline("Scaling background to aperture.")
    curve['bg_est'] = (np.array(curve['bg']) *
            (mc.area(radius) / (mc.area(annulus[1])-mc.area(annulus[0])))
            if annulus else np.array(curve['bg']))
    # FIXME: This exposure time calculation is a bottleneck!
    if verbose:
        mc.print_inline("Computing exposure times... could take a while...")
    curve['exptime'] = np.array(
            [dbt.compute_exptime(band, [curve['t0'][i], curve['t1'][i]])#, skypos=[ra0,dec0])
                for i in range(len(curve['t0']))])
    curve['exptime_data'] = np.array(
                [dbt.compute_exptime(band, [curve['t0_data'][i], curve['t1_data'][i]])#, skypos=[ra0,dec0])
                for i in range(len(curve['t0_data']))])
    curve['cps'] = (curve['counts']-curve['bg_est']) / curve['exptime']
    if verbose:
        mc.print_inline("Done.")
        mc.print_inline("")
    return curve

