"""Generate plots of difference between MCAT and gAperture photometry of
sources (delta-magnitude) as a function of source brightness (magnitude).

The following terminal commands will generate files containing gAperture and
MCAT photometry values for a large number of random-ish sources. They will
take a while to run (like a few days).

    ./gCalrun -f 'DPFCore_calrun_FUV.csv' -b 'FUV' --rarange [0,360] --decrange [-90,90] -n 150 --seed 323 -v 1

    ./gCalrun -f 'DPFCore_calrun_NUV.csv' -b 'NUV' --rarange [0,360] --decrange [-90,90] -n 150 --seed 323 -v 1
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plot
from matplotlib import gridspec

from gPhoton.regtestutils import datamaker
import gPhoton.galextools as gt
import gPhoton.MCUtils as mc
import gPhoton.gphoton_utils as gu
import gPhoton.gFind
import matplotlib.pyplot as plt

import astropy.coordinates as coord
import astropy.units as u

from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV

# Define I/O directories
outpath = '.'
print 'Writing to: {o}'.format(o=outpath)
inpath = '.'
print 'Reading from: {i}'.format(i=inpath)

# Setup reference values
scl = 1.4
bands = ['NUV','FUV']
base = 'DPFCore_calrun_'
apertxt='aper4'

# Read in the data
data = {}
for band in bands:
    filename = '{path}/{base}{band}.csv'.format(path=inpath,base=base,band=band)
    print filename
    data[band] = pd.read_csv(filename)#,nrows=None if band is 'FUV' else 10000)
    print '{band} sources: {cnt}'.format(
                                band=band,cnt=data[band]['objid'].shape[0])

# The following plot will demonstrate that there is good sky sampling
#for band in bands:
#    ra = coord.Angle(data[band]['ra']*u.degree)
#    ra = ra.wrap_at(180*u.degree)
#    dec = coord.Angle(data[band]['dec']*u.degree)
#    plt.title(band)
#    fig = plt.figure(figsize=(8,6))
#    plt.title(band)
#    ax = fig.add_subplot(111, projection="mollweide")
#    ax.scatter(ra.radian, dec.radian)

def make_kde(data,datarange,bwrange=[0.01,1]):
    # A function for producing Kernel Density Estimates
    # Based on code from:
    #   https://jakevdp.github.io/blog/2013/12/01/kernel-density-estimation/
    print 'Building KDE.'
    grid = GridSearchCV(KernelDensity(),
        {'bandwidth': np.linspace(bwrange[0],bwrange[1],100)},cv=20,n_jobs=4)
    grid.fit(data[:, None])
    bandwidth = grid.best_params_['bandwidth']
    x = np.linspace(datarange[0],datarange[1],10000)
    kde_skl = KernelDensity(bandwidth=bandwidth)
    kde_skl.fit(data[:, np.newaxis])
    # score_samples() returns the log-likelihood of the samples
    log_pdf = kde_skl.score_samples(x[:, np.newaxis])
    y = np.exp(log_pdf)
    peak = x[np.where(y==y.max())][0]
    return x,y,peak,bandwidth

# Generate dmag v. mag plots in both bands using both the Annulus and MCAT
# background estimation methods.
bincnt = 50
magstep = 1.
maglimits = [14,22.5]
magrange = np.arange(maglimits[0],maglimits[1]+magstep,magstep)
for band in bands:
    magmedian = np.zeros(len(magrange)-1)
    dmag = {#'NoBg':data[band]['aper4']-data[band]['mag'],
            'Annulus':data[band][apertxt]-data[band]['mag_bgsub'],
            'MCAT':data[band][apertxt]-data[band]['mag_mcatbgsub']}
    for bgmode in dmag.keys():
        dmagrange = [-0.5,0.5]
        fig = plt.figure(figsize=(8*scl,4*scl))
        fig.subplots_adjust(
            left=0.12,right=0.95,wspace=0.05,bottom=0.15,top=0.9)
        plt.subplot(1,2,1)
        ix = ((np.bitwise_and(np.array(data[band]['flags'].values,
                dtype='int16'),0b00111111)==0) & (data[band][apertxt]>0) &
                (data[band][apertxt]<=maglimits[1]))
        plt.xlabel('{b} AB Magnitude (MCAT {at})'.format(
                                b=band,at=str.upper(apertxt)),fontsize=14)
        plt.ylabel('{d}Magnitude (MCAT {at} - gAperture)'.format(
            d=r'$\Delta$',at=str.upper(apertxt)),fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.xlim([14,22.5])
        plt.ylim(dmagrange)
        # Compute moving median value in 1-mag bins
        for i,m in enumerate(magrange[:-1]):
            mix = ((np.bitwise_and(np.array(data[band]['flags'].values,
                dtype='int16'),0b00111111)==0) & (data[band][apertxt]>=m) &
                (data[band][apertxt]<m+magstep) & np.isfinite(dmag[bgmode]))
            magmedian[i]=dmag[bgmode][mix].median()
        plt.plot(data[band][apertxt][ix],dmag[bgmode][ix],'.',color='k',
            alpha=0.1 if band is 'FUV' else 0.05)
        plt.plot(magrange[:-1]+magstep/2.,magmedian,color='r',
                                            linestyle='dashed',linewidth=4)
        plt.subplot(1,2,2,ylim=dmagrange,yticks=[],xticks=[])
        plt.hist(np.array(dmag[bgmode][ix]),bins=bincnt,range=dmagrange,
                 orientation='horizontal',color='k',histtype='step',normed=1)
        x,y,peak,bandwidth = make_kde(
            np.array(dmag[bgmode][ix])[np.where(np.isfinite(dmag[bgmode][ix]))],
            dmagrange,bwrange=[0.01,0.1])
        print '{bgmode} ({b}): peak={p} +/- {p90} ({bw})'.format(
            bgmode=bgmode, b=band, p=peak,
            p90=np.percentile(np.abs(np.array(
                dmag[bgmode][ix]))[np.where(np.isfinite(dmag[bgmode][ix]))],90),
                bw=bandwidth)
        plt.plot(y,x) # Flipped to cheat the horizontal rotation
        plt.axhline(peak, color='k', linestyle='dotted', linewidth=2,
            label='KDE Peak: {p}'.format(p='{0:.2f}'.format(round(peak,2))))
        plt.axhline(dmag[bgmode][ix].median(),
            color='r', linestyle='dashed',linewidth=4,
            label='Median: {m}'.format(m=round(
                dmag[bgmode][ix].median(),2)))
        plt.text(0.6, -0.4, '50% of data within {pm}{p90}'.format(pm=r'$\pm$',
            p90=round(np.percentile(np.abs(np.array(
                dmag[bgmode][ix]))[np.where(
                                np.isfinite(dmag[bgmode][ix]))],50),2)),
            fontsize=15)
        plt.legend(fontsize=14)
        fig.savefig('{path}/Fig0{n}{l}.pdf'.format(path=outpath,
            n='4' if bgmode is 'Annulus' else '5',
            l='a' if band is 'NUV' else 'b',
            band=band,bg=bgmode.lower()),format='pdf',dpi=1000,
            bbox_inches='tight')

# Generate relative distributions (delta-mag) for background surface
# brightnesses estimated using the annulus vs. MCAT methods.
bincnt = 50
fig = plt.figure(figsize=(8*scl,4*scl))
fig.subplots_adjust(left=0.12,right=0.95,wspace=0.1,bottom=0.15,top=0.9)
for i,band in enumerate(bands):
    dmagrange = [-0.4,0.05]
    gphot_bg = data[band]['bg']/data[band]['exptime']
    mcat_bg = data[band]['mcat_bg.1']
    delta = mcat_bg - gphot_bg
    ix = (np.bitwise_and(np.array(data[band]['flags'].values,
            dtype='int16'),0b00111111)==0)
    plt.subplot(1,2,i+1,yticks=[])
    plt.hist(delta[ix],bins=bincnt,range=dmagrange,color='k',histtype='step',
        normed=1)
    x,y,peak,bandwidth = make_kde(delta[ix],dmagrange)
    plt.plot(x,y)
    plt.axvline(peak, color='k', linestyle='dotted', linewidth=2,
        label='KDE Peak: {p}'.format(p=round(peak,2)))
    plt.axvline(np.median(delta[ix]), color='r', linestyle='dashed',
        linewidth=4,
        label='Median: {m}'.format(m=round(np.median(delta[ix]),2)))
    plt.xlim(dmagrange)
    plt.xlabel('{b} {d}Magnitude/arcsec{exp} (MCAT - gAperture)'.format(
        b=band,d=r'$\Delta$',exp=r'$^{2}$'),fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.legend(loc=2,fontsize=12)
fig.savefig('{path}/Fig03.pdf'.format(path=outpath),
    format='pdf',dpi=1000,bbox_inches='tight')
#plt.close()

# Generate plots of relative astrometry between MCAT and gAperture using
# the center of brightness (CoB) from gAperture and reported source positions
# from the MCAT
bincnt = 101
a = 3600
dasrange = [-6,6]
for i,band in enumerate(bands):
    fig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[3,1])
    gs.update(left=0.1,right=0.95,wspace=0.02,hspace=0.02,bottom=0.1,top=0.95)
    delta_ra = np.array(
        data[band]['ra']-data[band]['racent'])*np.cos(np.array(data[band]['dec']))*a
    delta_dec = np.array(data[band]['dec']-data[band]['deccent'])*a
    ix = np.where((np.bitwise_and(np.array(data[band]['flags'].values,
            dtype='int16'),0b00111111)==0) & (data[band][apertxt]>0) &
                            (data[band][apertxt]<=maglimits[1]) &
                            np.isfinite(delta_ra) &
                            np.isfinite(delta_dec))
    plt.subplot(gs[2],xlim=dasrange,yticks=[])
    plt.gca().invert_yaxis()
    plt.xlabel('{d} Right Ascension (arcseconds)'.format(
        d=r'$\Delta$'),fontsize=14)
    plt.hist(delta_ra[ix],bins=bincnt,histtype='step',range=dasrange,
        normed=1,color='k')
    x,y,peak,bandwidth = make_kde(delta_ra[ix],dasrange)
    plt.plot(x,y)
    plt.axvline(peak, color='k', linestyle='dotted', linewidth=2,
        label='KDE Peak: {p} as'.format(p='{0:.2f}'.format(round(peak,2))))
    plt.axvline(np.median(delta_ra[ix]), color='r', linestyle='dashed',
        linewidth=4,
        label='Median: {m} as'.format(m=round(np.median(delta_ra[ix]),2)))
    ra_kdepeak, ra_median = peak, np.median(delta_ra[ix])
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.subplot(gs[1],ylim=dasrange,xticks=[],yticks=[])
    plt.hist(delta_dec[ix],bins=bincnt,histtype='step',
        orientation='horizontal',range=dasrange,normed=1,color='k')
    x,y,peak,bandwidth = make_kde(delta_dec[ix],dasrange)
    plt.plot(y,x)
    plt.axhline(peak, color='k', linestyle='dotted', linewidth=2,
        label='KDE Peak: {p} as'.format(p='{0:.2f}'.format(round(peak,2))))
    plt.axhline(np.median(delta_dec[ix]), color='r', linestyle='dashed',
        linewidth=4,
        label='Median: {m} as'.format(m=round(np.median(delta_dec[ix]),2)))
    dec_kdepeak, dec_median = peak, np.median(delta_dec[ix])
    plt.subplot(gs[0],xticks=[])
    plt.title('{b}: Relative Astrometry (MCAT - gAperture)'.format(
        b=band),fontsize=16)
    plt.plot(delta_ra[ix],delta_dec[ix],'.',color='k',alpha=0.1)
    plt.xlim(dasrange[0],dasrange[1])
    plt.ylim(dasrange[0],dasrange[1])
    plt.text(-5,4,'RA KDE Peak: {ra_kdepeak}\nRA Median: {ra_median}\nDec KDE Peak: {dec_kdepeak}\nDec Median: {dec_median}'.format(
            ra_kdepeak='{0:.2f}'.format(round(ra_kdepeak,3)),
            ra_median='{0:.2f}'.format(round(ra_median,3)),
            dec_kdepeak='{0:.2f}'.format(round(dec_kdepeak,3)),
            dec_median='{0:.2f}'.format(round(dec_median,))),
        bbox=dict(boxstyle='square',facecolor='w',edgecolor='k'))
    plt.ylabel('{d} Declination (arcseconds)'.format(d=r'$\Delta$'),fontsize=16)
    fig.savefig('{path}/Fig02{l}.pdf'.format(
        path=outpath,l='a' if band is 'NUV' else 'b'),
        format='pdf',dpi=1000,bbox_inches='tight')

#plt.close('all')
