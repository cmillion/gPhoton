"""Generate plots of difference between MCAT and gAperture photometry of
sources (delta-magnitude) as a function of source brightness (magnitude).

The following terminal commands will generate files containing gAperture and
MCAT photometry values for a large number of random-ish sources. They will
take a while to run (like a few days).

    ./gCalrun -f 'DPFCore_calrun_FUV.csv' -b 'FUV' --rarange [0,360] --decrange [-90,90] -n 150 --seed 323 -v 1

    ./gCalrun -f 'DPFCore_calrun_NUV.csv' -b 'NUV' --rarange [0,360] --decrange [-90,90] -n 150 --seed 323 -v 1
"""

import numpy as np
import matplotlib.pyplot as plot

from gPhoton.regtestutils import datamaker
import gPhoton.galextools as gt
import gPhoton.MCUtils as mc
import pandas as pd
import gPhoton.gphoton_utils as gu
import gPhoton.gFind
import numpy as np
import matplotlib.pyplot as plt

from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV

# Define I/O directories
outpath = '.'#'calpaper/paper'
inpath = '.'

# Setup reference values
scl = 1.4
bands = ['NUV','FUV']
base = 'DPFCore_calrun_'

# Read in the data
data = {}
for band in bands:
    filename = '{path}/{base}{band}.csv'.format(path=inpath,base=base,band=band)
    print filename
    data[band] = pd.read_csv(filename)
    print '{band} sources: {cnt}'.format(
                                band=band,cnt=data[band]['objid'].shape[0])

def make_kde(data,datarange,bwrange=[0.01,1]):
    # A function for producing Kernel Density Estimates
    # Based on code from:
    #   https://jakevdp.github.io/blog/2013/12/01/kernel-density-estimation/
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
magrange = np.arange(14,25,1)
for band in bands:
    magmedian = np.zeros(len(magrange)-1)
    dmag = {#'NoBg':data[band]['aper4']-data[band]['mag'],
            'Annulus':data[band]['aper4']-data[band]['mag_bgsub'],
            'MCAT':data[band]['aper4']-data[band]['mag_mcatbgsub']}
    for bgmode in dmag.keys():
        dmagrange = [-1,1]
        fig = plt.figure(figsize=(8*scl,4*scl))
        fig.subplots_adjust(
            left=0.12,right=0.95,wspace=0.05,bottom=0.15,top=0.9)
        plt.subplot(1,2,1)
        ix = ((data[band]['flags']==0) & (data[band]['aper4']>0))
        plt.title('{band}: {d}Mag vs. Mag ({bgmode} BG, n={n})'.format(
            band=band,d=r'$\Delta$',bgmode=bgmode,n=ix.shape[0]),fontsize=14)
        plt.xlabel('AB Magnitude (MCAT APER4)',fontsize=14)
        plt.ylabel('{d}Magnitude (MCAT APER4 - gAperture)'.format(
            d=r'$\Delta$'),fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.xlim([14,24])
        plt.ylim(dmagrange)
        for i,m in enumerate(magrange[:-1]):
            mix = ((data[band]['flags']==0) & (data[band]['aper4']>=m) &
                (data[band]['aper4']<m+1))
            magmedian[i]=dmag[bgmode][mix].median()
        plt.plot(data[band]['aper4'][ix],dmag[bgmode][ix],'.',color='k',
            alpha=0.1 if band is 'FUV' else 0.05)
        plt.plot(magrange[:-1]+0.5,magmedian,color='r',linestyle='dashed',
            linewidth=4)
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
        plt.axhline(np.median(dmag[bgmode][ix]), color='r', linestyle='dashed',
            linewidth=4,
            label='Median: {m}'.format(m=round(np.median(dmag[bgmode][ix]),2)))
        plt.text(0.5, -0.7, '50% of data within {pm}{p90}'.format(pm=r'$\pm$',
            p90=round(np.percentile(np.abs(np.array(
            dmag[bgmode][ix]))[np.where(np.isfinite(dmag[bgmode][ix]))],50),2)),
            fontsize=15)
        plt.legend(fontsize=14)
        fig.savefig('{path}/FigRelPhot{band}-{bg}_bg.pdf'.format(path=outpath,
            band=band,bg=bgmode.lower()),format='pdf',dpi=1000,
            bbox_inches='tight')
        plt.close()

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
    ix = (data[band]['flags']==0)
    plt.subplot(1,2,i+1,yticks=[])
    plt.title('{band} Relative BG (n={n})'.format(
        band=band,n=ix.shape[0]),fontsize=16)
    plt.hist(delta[ix],bins=bincnt,range=dmagrange,color='k',histtype='step',
        normed=1)
    x,y,peak,bandwidth = make_kde(delta[ix],dmagrange)
    print '{b}: peak={p} +/- {p90} ({bw})'.format(b=band,p=peak,bw=bandwidth,
        p90=np.percentile(np.abs(np.array(
        delta[ix]))[np.where(np.isfinite(dmag[bgmode][ix]))],90))
    plt.plot(x,y)
    plt.axvline(peak, color='k', linestyle='dotted', linewidth=2,
        label='KDE Peak: {p}'.format(p=round(peak,2)))
    plt.axvline(np.median(delta[ix]), color='r', linestyle='dashed',
        linewidth=4,
        label='Median: {m}'.format(m=round(np.median(delta[ix]),2)))
    plt.xlim(dmagrange)
    plt.xlabel('{d}Magnitude/arcsec{exp} (MCAT - gAperture)'.format(
        d=r'$\Delta$',exp=r'$^{2}$'),fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.legend(loc=2,fontsize=12)
fig.savefig('{path}/BG_dMag.pdf'.format(path=outpath),
    format='pdf',dpi=1000,bbox_inches='tight')
plt.close()

# Generate plots of relative astrometry between MCAT and gAperture using
# the center of brightness (CoB) from gAperture and reported source positions
# from the MCAT
bincnt = 50
a = 3600
dasrange = [-6,6]
for i,band in enumerate(bands):
    fig = plt.figure(figsize=(8*scl,8*scl))
    fig.subplots_adjust(left=0.12,right=0.95,wspace=0.02,hspace=0.02,
        bottom=0.15,top=0.9)
    delta_ra = np.array(
        data[band]['ra']-data[band]['racent'])*np.cos(np.array(data[band]['dec']))*a
    delta_dec = np.array(data[band]['dec']-data[band]['deccent'])*a
    ix = np.where((data[band]['flags']==0) & (data[band]['ra']>1) &
        (data[band]['ra']<359) & np.isfinite(delta_ra) & np.isfinite(delta_dec))
    plt.subplot(2,2,1,xticks=[])
    plt.title('{b}: Relative Astrometry (MCAT - gAperture)'.format(
        b=band),fontsize=16)
    plt.plot(delta_ra[ix],delta_dec[ix],'.',color='k',alpha=0.1)
    plt.xlim(dasrange[0],dasrange[1])
    plt.ylim(dasrange[0],dasrange[1])
    plt.ylabel('{d} Declination (arcseconds)'.format(d=r'$\Delta$'),fontsize=14)
    plt.subplot(2,2,3,xlim=dasrange,yticks=[])
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
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.legend(loc=3,fontsize=14)
    plt.subplot(2,2,2,ylim=dasrange,xticks=[],yticks=[])
    plt.hist(delta_dec[ix],bins=bincnt,histtype='step',
        orientation='horizontal',range=dasrange,normed=1,color='k')
    x,y,peak,bandwidth = make_kde(delta_dec[ix],dasrange)
    plt.plot(y,x)
    plt.axhline(peak, color='k', linestyle='dotted', linewidth=2,
        label='KDE Peak: {p} as'.format(p='{0:.2f}'.format(round(peak,2))))
    plt.axhline(np.median(delta_dec[ix]), color='r', linestyle='dashed',
        linewidth=4,
        label='Median: {m} as'.format(m=round(np.median(delta_dec[ix]),2)))
    plt.legend(fontsize=14)
    fig.savefig('{path}/RelAstrometry{band}.pdf'.format(path=outpath,band=band),
        format='pdf',dpi=1000,bbox_inches='tight')
    plt.close()
