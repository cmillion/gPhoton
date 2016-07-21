"""Generate plots of difference between MCAT and gAperture photometry of
sources (delta-magnitude) as a function of source brightness (magnitude).

The following terminal commands will generate files containing gAperture and
MCAT photometry values for a large number of random-ish sources. They will
take a while to run (like a few days).

    ./gCalrun -f 'DPFCore_calrun_FUV.csv' -b 'FUV' --rarange [0,360] --decrange [-90,90] -n 40000 --seed 323 -v 1

    ./gCalrun -f 'DPFCore_calrun_NUV.csv' -b 'NUV' --rarange [0,360] --decrange [-90,90] -n 40000 --seed 323 -v 1

The number of 'draws' defined by `-n` has been set to an arbitrarily high
value. Only the first 10,000 entries in the output files are used in the
subsequent analyses defined by the scripts that follow.
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
    data[band] = pd.read_csv(filename,nrows=10000)
    print '{band} sources: {cnt}'.format(
                                band=band,cnt=data[band]['objid'].shape[0])

# The following plot will demonstrate that there is good sky sampling
for band in bands:
   ra = coord.Angle(data[band]['ra']*u.degree)
   ra = ra.wrap_at(180*u.degree)
   dec = coord.Angle(data[band]['dec']*u.degree)
   plt.title(band)
   fig = plt.figure(figsize=(8,6))
   plt.title(band)
   ax = fig.add_subplot(111, projection="mollweide")
   ax.scatter(ra.radian, dec.radian)

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
print 'dMag v. Mag'
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
        print 'Using {n} of {m} {b} sources.'.format(n=np.where(ix)[0].size,
            m=np.array(data[band]['flags']).size,b=band)
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
            band=band,bg=bgmode.lower()),format='pdf',dpi=1000)#,bbox_inches='tight')

# Generate relative distributions (delta-mag) for background surface
# brightnesses estimated using the annulus vs. MCAT methods.
bincnt = 50
fig = plt.figure(figsize=(8*scl,4*scl))
fig.subplots_adjust(left=0.12,right=0.95,wspace=0.1,bottom=0.15,top=0.9)
print 'bg dMag v. Mag (surface)'
for i,band in enumerate(bands):
    dmagrange = [-0.4,0.05]
    gphot_bg = data[band]['bg']/data[band]['exptime']
    mcat_bg = data[band]['mcat_bg.1']
    delta = mcat_bg - gphot_bg
    ix = (np.bitwise_and(np.array(data[band]['flags'].values,
            dtype='int16'),0b00111111)==0)
    print 'Using {n} of {m} {b} sources.'.format(n=np.where(ix)[0].size,
            m=np.array(data[band]['flags']).size,b=band)
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
    format='pdf',dpi=1000)#,bbox_inches='tight')
#plt.close()

# Generate plots of relative astrometry between MCAT and gAperture using
# the center of brightness (CoB) from gAperture and reported source positions
# from the MCAT
bincnt = 101
a = 3600
dasrange = [-6,6]
print 'Astrometry'
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
    print 'Using {n} of {m} {b} sources.'.format(n=np.where(ix)[0].size,
            m=np.array(data[band]['flags']).size,b=band)
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
        format='pdf',dpi=1000)#,bbox_inches='tight')

#plt.close('all')

#############################

""" Generate absolute photometry plots of the GALEX white dwarf standard star
LDS749B.

The following commands should be copied and pasted into a Python terminal w/
gPhoton installed in order to generate the source data. The function thus
called queries the MCAT for all GALEX observations of LDS749B and constructs
gAperture calls equivalent sky positions and time ranges. The MCAT and
gAperture data are then compiled into a single reference file for comparison.
There is a lot of data, so these will take a long while to run. The functions
write as they go and check objids against what has already been processed, so
you can resume crashed or canceled jobs as long as the output file is the same.

LDS749B w/ 34.5 arcsecond aperture & large annulus

    from gPhoton.regtestutils import datamaker
    datamaker('FUV',[323.06766667,0.25400000],'LDS749B_dm_FUV_35as.csv',margin=0.001,searchradius=0.001,annulus=[0.025,0.05],radius=0.004805555555555556)

    from gPhoton.regtestutils import datamaker
    datamaker('NUV',[323.06766667,0.25400000],'LDS749B_dm_NUV_35as.csv',margin=0.001,searchradius=0.001,annulus=[0.025,0.05],radius=0.004805555555555556)

The following commands will generate the plots. In addition to gPhoton and its
dependencies, you will need matplotlib and scikit-learn.
"""
import numpy as np
import matplotlib.pyplot as plt

import gPhoton.galextools as gt
import gPhoton.dbasetools as db
import gPhoton.gQuery as gq
import gPhoton.MCUtils as mc
from gPhoton.gphoton_utils import read_lc, dmag_errors, data_errors

from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV

# Define I/O directories
outpath = '.'
print 'Writing to: {o}'.format(o=outpath)
inpath = '.'
print 'Reading from: {i}'.format(i=inpath)

# Setup reference variables.
scl = 1.4
bands = ['FUV','NUV']

# Read in the calibration data
source = 'LDS749B'
aperas = '35as'
target = '{s}_{a}'.format(s=source,a=aperas)
refmag = {'FUV':15.6,'NUV':14.76} # per http://arxiv.org/pdf/1312.3882.pdf
#refmag = {'FUV':15.57,'NUV':14.71} # per calibration paper
#refmag = {'FUV':15.67,'NUV':14.82} # GALEX "published" (?) per Camarota, et al.
data = {'FUV':read_lc('{path}/{s}_dm_FUV_{a}.csv'.format(
                                            s=source,path=inpath,a=aperas)),
        'NUV':read_lc('{path}/{s}_dm_NUV_{a}.csv'.format(
                                            s=source,path=inpath,a=aperas))}
aper = 7
radius = gt.aper2deg(aper)
apertxt = 'aper{n}'.format(n=int(aper))

# To make a cut on AIS leg, run the following and then add the leg
# condition in ix below.
try:
    _ = len(data['FUV']['legs'])
except KeyError:
    print 'Retrieving observation leg information from database.'
    data['FUV']['legs'] = np.array([gq.getArray(
            gq.obstype(o))[0][5] for o in np.array(data['FUV']['objid'])])
    data['FUV'].to_csv('{path}/{s}_dm_FUV_{a}.csv'.format(
                                                s=source,path=inpath,a=aperas))
legs = {'FUV':np.array(data['FUV']['legs'])}

# Provide some baseline statistics for the sample
ix = {}
for band in bands:
    # The bitwise_and ignores the annulus mask flags because using MCAT bg.
    ix[band] = np.where((data[band][apertxt]>0) &
                        (data[band]['detrad']<200) &
                        (data[band]['t0']<961986575.) &
        (np.bitwise_and(np.array(data[band]['flags'].values,dtype='int16'),
                                        0b00111111)==0))# & (legs[band]>3)
    print '{b}: {m} / {n} data points used'.format(b=band,
        n=len(data[band][apertxt]),m=len(ix[band][0]))
    print 'mcat: {m} (ref: {r})'.format(
        m = np.median(np.array(data[band][apertxt])[ix[band]]) -
                                                    gt.apcorrect1(radius,band),
        r = refmag[band])

# The following plots will demonstrate the FUV multi-modality
#for band in bands:
fig = plt.figure(figsize=(10*scl,4*scl))
fig.subplots_adjust(left=0.12,right=0.95,wspace=0.05,bottom=0.15,top=0.9)
band = 'FUV'
plt.subplot(1,2,1)
#plt.title('{s}: {b} Multi-modality'.format(s=source,b=band))
plt.ylabel('{d} AB Magnitude (MCAT - gAperture)'.format(d=r'$\Delta$'))
plt.xlabel('CAI Observation Leg ({b}, {target})'.format(
                                    target=target.split('_')[0],b=band))
plt.axhline(-0.015, color='g', linestyle='solid', linewidth=1)
plt.axhline(-0.05, color='g', linestyle='solid', linewidth=1)
dmag = data[band]['aper{a}'.format(a=aper)]-data[band]['mag_mcatbgsub']
dmagrange = [-0.2,0.025]
plt.ylim(dmagrange)
plt.plot(legs[band][ix[band]],
                        np.array(dmag)[ix[band]],'.',alpha=0.3,color='k')
plt.subplot(1,2,2,yticks=[])
plt.hist(np.array(dmag)[ix[band]],bins=50,range=dmagrange,
                 orientation='horizontal',color='k',histtype='step',normed=1)
plt.axhline(-0.015, color='g', linestyle='solid', linewidth=1)
plt.axhline(-0.05, color='g', linestyle='solid', linewidth=1)
plt.savefig('{path}/Fig09a.pdf'.format(path=outpath),
    format='pdf',dpi=1000)#,bbox_inches='tight')

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

# Overplot gAperture photometry of LDS749B on the reference magnitude and
# modeled 3-sigma error bounds as a function of exposure time.
nsigma = 3
for band in data.keys():
    tmin,tmax = 0,300
    plt.figure(figsize=(8*scl,4*scl))
    plt.subplots_adjust(left=0.12,right=0.95,wspace=0.1,bottom=0.15,top=0.9)
    # Overplot data in AIS legs 1-3.
    if band is 'FUV':
        lix = np.where(legs[band][ix[band]]<=3)
        plt.errorbar(np.array(data[band]['exptime'])[ix[band]][lix],
            np.array(data[band]['mag_mcatbgsub'])[ix[band][0][lix]]-gt.apcorrect1(radius,band),
            fmt='x',color='r',alpha=0.3,
            yerr=[np.array(data[band]['mag_mcatbgsub_err_1'])[ix[band][0][lix]],
                np.array(data[band]['mag_mcatbgsub_err_2'])[ix[band][0][lix]]])
    plt.errorbar(np.array(data[band]['exptime'])[ix[band]],
        np.array(data[band]['mag_mcatbgsub'])[ix[band]]-gt.apcorrect1(radius,band),
        fmt='.',color='k',alpha=0.1,
        yerr=[np.array(data[band]['mag_mcatbgsub_err_1'])[ix[band]],
              np.array(data[band]['mag_mcatbgsub_err_2'])[ix[band]]])
    t = np.arange(tmin+1,tmax+1)
    plt.plot(t,refmag[band]*(t/t),'k')
    plt.plot(t,data_errors(refmag[band],band,t,sigma=nsigma)[0],'r--')
    plt.plot(t,data_errors(refmag[band],band,t,sigma=nsigma)[1],'r--')
    plt.xlim(tmin,tmax)
    plt.ylim(refmag[band]-0.5,refmag[band]+1)
    plt.gca().invert_yaxis()
    plt.xlabel('Effective Exposure Depth (s, n={n})'.format(
        n=len(ix[band][0])),fontsize=16)
    plt.ylabel('{b} gAperture Magnitude ({target})'.format(
        target=target.split('_')[0],b=band),fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    a,b = data_errors(refmag[band],band,
        np.array(data[band]['exptime'])[ix[band]],sigma=nsigma)
    cnt = len(np.where(
        (np.array(data[band]['mag_mcatbgsub'])[ix[band]]-
         gt.apcorrect1(radius,band)>=b) &
        (np.array(data[band]['mag_mcatbgsub'])[ix[band]]-
         gt.apcorrect1(radius,band)<=a))[0])
    print '{b}: {n} of {m} ({p}%) within {s} sigma'.format(
        b=band,n=cnt,m=len(ix[band][0]),p=100*cnt/len(ix[band][0]),s=nsigma)
    if band is 'FUV':
        lix = np.where((legs[band][ix[band]]>3) | (legs[band][ix[band]]<0))
        cntlix = len(np.where(
            (np.array(data[band]['mag_mcatbgsub'])[ix[band][0][lix]]-
                 gt.apcorrect1(radius,band)>=b[lix]) &
            (np.array(data[band]['mag_mcatbgsub'])[ix[band][0][lix]]-
                 gt.apcorrect1(radius,band)<=a[lix]))[0])
        print 'w/o legs 1-3... {b}: {n} of {m} ({p}%) within {s} sigma'.format(
            b=band,n=cntlix,m=len(ix[band][0][lix]),p=100*cntlix/len(ix[band][0][lix]),s=nsigma)
    plt.text(150, 16.3 if band=='FUV' else 15.4, '{p}% within {s}{sym}'.format(
        p=100*cnt/len(ix[band][0]),s=nsigma,sym=r'$\sigma$'), fontsize=18)
    plt.savefig('{path}/Fig0{n}a.pdf'.format(path=outpath,
                n='6' if band is 'NUV' else '7'),
                format='pdf',dpi=1000)#,bbox_inches='tight')

# Generate distribution plots for gAperture photometry of LDS749B
bincnt=50
fig = plt.figure(figsize=(10*scl,4*scl))
fig.subplots_adjust(left=0.12,right=0.95,wspace=0.1,bottom=0.15,top=0.9)
for i,band in enumerate(['NUV','FUV']):
    magrange = [refmag[band]-0.2,refmag[band]+0.7]
    tmin,tmax = 0,350
    plt.subplot(1,2,i+1,xticks=np.arange(round(magrange[0],1),
        round(magrange[1]+0.01,1),0.2),yticks=[])
    mags = np.array(data[band]['mag_mcatbgsub'])-gt.apcorrect1(radius,band)
    plt.hist(mags[ix[band]],bins=bincnt,color='k',histtype='step',
        range=magrange,normed=1)
    plt.axvline(refmag[band], color='g', linestyle='solid', linewidth=4,
        label='Ref: {r} AB Mag'.format(r=round(refmag[band],2)))
    x,y,peak,bandwidth = make_kde(mags[ix[band]],magrange)
    plt.plot(x,y)
    print '{b}: peak={p} ({bw})'.format(b=band,p=peak,bw=bandwidth)
    if band is 'FUV':
        lix = np.where((legs[band][ix[band]]>3) | (legs[band][ix[band]]<0))
        _,_,newpeak,_ = make_kde(mags[ix[band][0][lix]],magrange)
        print 'w/o legs 1-3, FUV peak at {p}'.format(p=newpeak)
    plt.axvline(peak, color='k', linestyle='dotted', linewidth=2,
        label='KDE Peak: {p}'.format(p=round(peak,2)))
    plt.axvline(np.median(mags[ix[band]]), color='r', linestyle='dashed',
        linewidth=4,
        label='Median: {m}'.format(m=round(np.median(mags[ix[band]]),2)))
    plt.xlim(magrange)
    plt.gca().invert_xaxis()
    plt.legend(loc=2,fontsize=14)
    plt.xlabel('{b} gAperture Magnitude ({target}, n={n})'.format(
        target=target.split('_')[0],b=band,n=len(ix[band][0])),fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.text(magrange[1]-0.1,3 if band is 'NUV' else 1,
                            '(a)' if band is 'NUV' else '(b)',fontsize=30)
    plt.savefig('{path}/Fig08a.pdf'.format(path=outpath),
        format='pdf',dpi=1000)#,bbox_inches='tight')

# Overplot MCAT photometry of LDS749B on the reference magnitude and
# modeled 3-sigma error bounds as a function of exposure time.
for band in data.keys():
    tmin,tmax = 0,300#data[band]['t_eff'].min(),data[band]['t_eff'].max()
    plt.figure(figsize=(8*scl,4*scl))
    plt.subplots_adjust(left=0.12,right=0.95,wspace=0.1,bottom=0.15,top=0.9)
    plt.errorbar(np.array(data[band]['exptime'])[ix[band]],
        np.array(data[band][apertxt])[ix[band]]-gt.apcorrect1(gt.aper2deg(aper),band),
            fmt='.',color='k',alpha=0.1,
            yerr=[np.array(data[band]['mag_bgsub_err_1'])[ix[band]],
                  np.array(data[band]['mag_bgsub_err_2'])[ix[band]]])
    t = np.arange(tmin+1,tmax+1)
    plt.plot(t,refmag[band]*(t/t),'k')
    plt.plot(t,data_errors(refmag[band],band,t,sigma=3)[0],'r--')
    plt.plot(t,data_errors(refmag[band],band,t,sigma=3)[1],'r--')
    plt.xlim(tmin,tmax)
    plt.ylim(refmag[band]-0.5,refmag[band]+1)
    plt.gca().invert_yaxis()
    plt.xlabel('Effective Exposure Depth (s, n={n})'.format(
        n=len(ix[band][0])),fontsize=16)
    plt.ylabel('{b} MCAT Magnitude ({target})'.format(
        target=target.split('_')[0],b=band),fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    a,b = data_errors(refmag[band],
        band,np.array(data[band]['exptime'])[ix[band]],sigma=nsigma)
    cnt = len(np.where((np.array(data[band][apertxt])[ix[band]]-
                                    gt.apcorrect1(gt.aper2deg(aper),band)>=b) &
                       (np.array(data[band][apertxt])[ix[band]]-
                                gt.apcorrect1(gt.aper2deg(aper),band)<=a))[0])
    print '{b}: {n} of {m} ({p}%) within {s} sigma'.format(
        b=band,n=cnt,m=len(ix[band][0]),p=100*cnt/len(ix[band][0]),s=nsigma)
    plt.text(150, 16.3 if band=='FUV' else 15.4, '{p}% within {s}{sym}'.format(
        p=100*cnt/len(ix[band][0]),s=nsigma,sym=r'$\sigma$'), fontsize=18)
    plt.savefig('{path}/Fig0{n}b.pdf'.format(path=outpath,
                n='6' if band is 'NUV' else '7'),
                format='pdf',dpi=1000)#,bbox_inches='tight')

# Generate distribution plots for MCAT photometry of LDS749B
bincnt = 50
fig = plt.figure(figsize=(10*scl,4*scl))
fig.subplots_adjust(left=0.12,right=0.95,wspace=0.1,bottom=0.15,top=0.9)
for i,band in enumerate(['NUV','FUV']):
    magrange = [refmag[band]-0.2,refmag[band]+0.6]
    tmin,tmax = 0,350
    plt.subplot(1,2,i+1,xticks=np.arange(round(magrange[0],1),
        round(magrange[1]+0.01,1),0.2),yticks=[])
    mags = np.array(data[band][apertxt])-gt.apcorrect1(gt.aper2deg(aper),band)
    plt.hist(mags[ix[band]],bins=bincnt,color='k',histtype='step',
        range=magrange, normed=1)
    plt.axvline(refmag[band], color='g', linestyle='solid', linewidth=4,
        label='Ref: {r} AB Mag'.format(r=round(refmag[band],2)))
    x,y,peak,bandwidth = make_kde(mags[ix[band]],magrange)
    print '{b}: peak={p} ({bw})'.format(b=band,p=peak,bw=bandwidth)
    plt.plot(x,y)
    plt.axvline(peak, color='k', linestyle='dotted', linewidth=2,
        label='KDE Peak: {p}'.format(p=round(peak,2)))
    plt.axvline(np.median(mags[ix[band]]), color='r', linestyle='dashed',
        linewidth=4,
        label='Median: {m}'.format(m=round(np.median(mags[ix[band]]),2)))
    plt.xlim(magrange)
    plt.gca().invert_xaxis()
    plt.xlabel('{b} MCAT Magnitude ({target}, n={n})'.format(
        target=target.split('_')[0],b=band,n=len(ix[band][0])),fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.legend(loc=2,fontsize=14)
    plt.text(magrange[1]-0.1,3 if band is 'NUV' else 2,
                            '(c)' if band is 'NUV' else '(d)',fontsize=30)
    plt.savefig('{path}/Fig08b.pdf'.format(path=outpath),
        format='pdf',dpi=1000)#,bbox_inches='tight')

####

# Overplot gAperture photometry of LDS749B on the reference magnitude and
# modeled 3-sigma error bounds as a function of exposure time
# EXCLUDING LEGS 1,2,3.
nsigma = 3
for band in ['FUV']:
    tmin,tmax = 0,300
    plt.figure(figsize=(8*scl,4*scl))
    plt.subplots_adjust(left=0.12,right=0.95,wspace=0.1,bottom=0.15,top=0.9)
    # Exclude data in AIS legs 1-3.
    lix = np.where(legs[band][ix[band]]>3)
    plt.errorbar(np.array(data[band]['exptime'])[ix[band]][lix],
        np.array(data[band]['mag_mcatbgsub'])[ix[band][0][lix]]-gt.apcorrect1(radius,band),
        fmt='.',color='k',alpha=0.1,
        yerr=[np.array(data[band]['mag_mcatbgsub_err_1'])[ix[band][0][lix]],
                np.array(data[band]['mag_mcatbgsub_err_2'])[ix[band][0][lix]]])
    t = np.arange(tmin+1,tmax+1)
    plt.plot(t,refmag[band]*(t/t),'k')
    plt.plot(t,data_errors(refmag[band],band,t,sigma=nsigma)[0],'r--')
    plt.plot(t,data_errors(refmag[band],band,t,sigma=nsigma)[1],'r--')
    plt.xlim(tmin,tmax)
    plt.ylim(refmag[band]-0.5,refmag[band]+1)
    plt.gca().invert_yaxis()
    plt.xlabel('Effective Exposure Depth (s, n={n})'.format(
        n=len(ix[band][0])),fontsize=16)
    plt.ylabel('{b} gAperture Magnitude ({target})'.format(
        target=target.split('_')[0],b=band),fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    a,b = data_errors(refmag[band],band,
        np.array(data[band]['exptime'])[ix[band]],sigma=nsigma)
    lix = np.where((legs[band][ix[band]]>3) | (legs[band][ix[band]]<0))
    cnt = len(np.where(
        (np.array(data[band]['mag_mcatbgsub'])[ix[band][0][lix]]-
             gt.apcorrect1(radius,band)>=b[lix]) &
        (np.array(data[band]['mag_mcatbgsub'])[ix[band][0][lix]]-
             gt.apcorrect1(radius,band)<=a[lix]))[0])
    print '{b}: {n} of {m} ({p}%) within {s} sigma'.format(
        b=band,n=cnt,m=len(ix[band][0]),p=100*cnt/len(ix[band][0]),s=nsigma)
    print 'w/o legs 1-3... {b}: {n} of {m} ({p}%) within {s} sigma'.format(
            b=band,n=cntlix,m=len(ix[band][0][lix]),p=100*cntlix/len(ix[band][0][lix]),s=nsigma)
    plt.text(150, 16.3 if band=='FUV' else 15.4, '{p}% within {s}{sym} (excluding legs 1-3)'.format(
        p=100*cntlix/len(ix[band][0][lix]),s=nsigma,sym=r'$\sigma$'), fontsize=18)
    plt.savefig('{path}/Fig09b.pdf'.format(path=outpath),
        format='pdf',dpi=1000)#,bbox_inches='tight')

# Generate distribution plots for gAperture photometry of LDS749B
# EXLUCIND LEGS 1,2,3
bincnt=50
fig = plt.figure(figsize=(5*scl,4*scl))
fig.subplots_adjust(left=0.12,right=0.95,wspace=0.1,bottom=0.15,top=0.9)
for i,band in enumerate(['FUV']):
    magrange = [refmag[band]-0.2,refmag[band]+0.7]
    tmin,tmax = 0,350
    lix = np.where(legs[band][ix[band]]>3)
    plt.subplot(1,1,i+1,xticks=np.arange(round(magrange[0],1),
        round(magrange[1]+0.01,1),0.2),yticks=[])
    mags = np.array(data[band]['mag_mcatbgsub'])-gt.apcorrect1(radius,band)
    plt.hist(mags[ix[band]][lix],bins=bincnt,color='k',histtype='step',
        range=magrange,normed=1)
    plt.axvline(refmag[band], color='g', linestyle='solid', linewidth=4,
        label='Ref: {r} AB Mag'.format(r=round(refmag[band],2)))
    x,y,peak,bandwidth = make_kde(mags[ix[band]][lix],magrange)
    plt.plot(x,y)
    print '{b}: peak={p} ({bw})'.format(b=band,p=peak,bw=bandwidth)
    if band is 'FUV':
        lix = np.where((legs[band][ix[band]]>3) | (legs[band][ix[band]]<0))
        _,_,newpeak,_ = make_kde(mags[ix[band][0][lix]],magrange)
        print 'w/o legs 1-3, FUV peak at {p}'.format(p=newpeak)
    plt.axvline(peak, color='k', linestyle='dotted', linewidth=2,
        label='KDE Peak: {p}'.format(p=round(peak,2)))
    plt.axvline(np.median(mags[ix[band]][lix]), color='r', linestyle='dashed',
        linewidth=4,
        label='Median: {m}'.format(m=round(np.median(mags[ix[band]][lix]),2)))
    plt.xlim(magrange)
    plt.gca().invert_xaxis()
    plt.legend(loc=2,fontsize=14)
    plt.text(16.2,3,'Excludes legs 1-3.',fontsize=16)
    plt.xlabel('{b} gAperture Magnitude ({target}, n={n})'.format(
        target=target.split('_')[0],b=band,n=len(ix[band][0])),fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.savefig('{path}/Fig09c.pdf'.format(path=outpath),
        format='pdf',dpi=1000)#,bbox_inches='tight')

#plt.close('all')

#########################
'''Stim rate / deadtime analyses.'''

"""Read and repackage the stim data that has been generated within the database
for a large number of eclipses and then generate linear fits to data in both
bands using MCMC in order to determine the stim rate as a function of global
count rate, which is a proxy for detector dead time.

The source file `stimquery.csv` was generated by running SQL queries _like_
the following for every entry in `FUVStim.csv` and `NUVStim.csv`. It was run on the server on account of being much, much faster.

    if object_id('tempdb.dbo.#temp1','U') IS NOT NULL drop table #temp1
    go
    select x, y into #temp1 from GPFCore.dbo.FUVPhotonsV
    where time >= 754607026995 and time < 754607158995
    union all
    select x, y from GPFCore.dbo.FUVPhotonsNULLV
    where time >= 754607026995 and time < 754607158995
    go
    create clustered index temp$1 on #temp1(x,y)
    go
    insert into stimResults select 2010,'150.910025185_75.9819518125','FUV',count(*) as theCount, getdate() as theTime from #temp1
    where (((x >= -38279.290277 and x< -35661.0053807) and (y >=34401.4593746 and y < 37019.7442709))
    or ((x >= 36972.6203787 and x< 39590.905275) and (y >=34397.6778277 and y < 37015.962724))
    or ((x >= -38274.199733 and x< -35655.9148368) and (y >=-38410.3354215 and y < -35792.0505252))
    or ((x >= 36967.0935025 and x< 39585.3783987) and (y >=-38410.7717538 and y < -35792.4868575)))
    go

The `stimquery.csv` file is read in and reformated as python pickle files
(.pkl). This extra and seemingly unnecessary step is a holdover from early
development. There is a lot of redundant code between the FUV and NUV mixture models, which is also just a holdover from early development.
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gPhoton.dbasetools as dt
import gPhoton.gQuery as gq
import pprint
import triangle
import emcee
import pickle

inpath = '.'
outpath = '.'

#------------

"""Read and repackage the stim data that has been generated within the database for a large number of eclipses."""
bands = ['FUV','NUV']
querydata = pd.read_csv('{path}/stimquery.csv'.format(path=inpath))

dtdata = {}
for band in bands:
    stimdata = pd.read_csv('{b}Stim.csv'.format(b=band),index_col=0)
    if len(stimdata.keys())==10:
        continue # The additional data has already been added to the file
    n = len(stimdata)
    dtdata[band] = {'t0':np.zeros(n),'t1':np.zeros(n),
                    'stimcount':np.zeros(n),'globalcount':np.zeros(n),
                    'exptime':np.zeros(n), 'deadtime':np.zeros(n)}
    for i,data in enumerate(np.array(stimdata)):
        ix = np.where(querydata['ix']==data[0])
        t0 = np.array(querydata['t0'][ix[0]])[0]/1000.
        t1 = np.array(querydata['t1'][ix[0]])[0]/1000.
        print i,data[0],band,t0,t1
        exptime = t1-t0-dt.compute_shutter(band,[t0,t1])
        stimcount = data[3]
        globalcount = gq.getValue(gq.globalcounts(band,t0,t1))
        deadtime = gq.getValue(gq.deadtime(band,t0,t1))
        dtdata[band]['t0'][i] = t0
        dtdata[band]['t1'][i] = t1
        dtdata[band]['stimcount'][i] = stimcount
        dtdata[band]['globalcount'][i] = globalcount
        dtdata[band]['exptime'][i] = exptime
        dtdata[band]['deadtime'][i] = deadtime

    for k in dtdata[band].keys():
        stimdata[k] = dtdata[band][k]
    stimdata.to_csv('{b}Stim.csv'.format(b=band))
#------------------------------------
# FUV linear mixture model

scl = 1.4

band = 'FUV'
dtdata = {'FUV':pd.read_csv(
    '{path}/FUVStim.csv'.format(path=inpath),index_col=0)}
print dtdata[band].keys()

ix = np.where(np.isfinite(dtdata[band]['globalcount']) &
    np.isfinite(dtdata[band]['stimcount']) & (dtdata[band]['exptime']>0) &
    (dtdata[band]['stimcount']/dtdata[band]['exptime']<140))
rawexpt = np.array(dtdata[band]['t1'])[ix]-np.array(dtdata[band]['t0'])[ix]
x = np.array(dtdata[band]['globalcount'])[ix]/rawexpt
y = np.array(dtdata[band]['stimcount'])[ix]/rawexpt
xerr = (np.sqrt(np.array(dtdata[band]['globalcount'])[ix])/
                                        np.array(dtdata[band]['exptime'])[ix])
yerr = (np.sqrt(np.array(dtdata[band]['stimcount'])[ix])/
                                        np.array(dtdata[band]['exptime'])[ix])
print x.min(),x.max(),xerr.min(),xerr.max()
print y.min(),y.max(),yerr.min(),yerr.max()

# Plots the raw data. Useful for diagnostics.
# plt.figure()
# plt.title('{b} Stim vs. Global Count Rates'.format(b=band))
# plt.xlabel('Global Count Rate (ct/s)')
# plt.ylabel('Stim Count Rate (ct/s)')
# plt.plot(x,y,'.',alpha=0.25,label='{b}'.format(b=band))
#
# plt.figure()
# plt.title('{b} Stim vs. Global Count Rates Errors'.format(b=band))
# plt.xlabel('Global Count Rate Error (ct/s)')
# plt.ylabel('Stim Count Rate Error (ct/s)')
# plt.plot(xerr,yerr,'.',alpha=0.25,label='{b}'.format(b=band))
#
# plt.show()

# Build a linear mixture model. Run MCMC on it.
# Based on http://dan.iel.fm/emcee/current/user/line/
# Define the probabilistic model...
def lnprior(p,bounds):
    # We'll just put reasonable uniform priors on all the parameters.
    if not all([b[0] < v < b[1] for v, b in zip(p, bounds)]):
        return -np.inf
    return 0

# The "foreground" linear likelihood:
def lnlike_fg(p,x,y,yerr):
    m, b, _, M, lnV = p
    model = m * x + b
    return -0.5 * (((model - y) / yerr) ** 2 + 2 * np.log(yerr))

# The "background" outlier likelihood:
def lnlike_bg(p,x,y,yerr):
    _, _, Q, M, lnV = p
    var = np.exp(lnV) + yerr**2
    return -0.5 * ((M - y) ** 2 / var + np.log(var))

# Full probabilistic model.
def lnprob(p,x,y,yerr,bounds):
    m, b, Q, M, lnV = p
    lp = lnprior(p,bounds)
    if not np.isfinite(lp):
        return -np.inf, None
    ll_fg = lnlike_fg(p,x,y,yerr)
    arg1 = ll_fg + np.log(Q)
    ll_bg = lnlike_bg(p,x,y,yerr)
    arg2 = ll_bg + np.log(1.0 - Q)
    ll = np.sum(np.logaddexp(arg1, arg2))
    return lp + ll, (arg1, arg2)

# Initialize the walkers at a reasonable location.
ndim, nwalkers = 5, 32
params = ['m', 'b', 'Q', 'M', 'lnV']
bounds = [(-.001,0), (75, 81), (0, 1), (0, 100), (0, 10)]
p0 = np.array([-0.0005, 79, 0.5, 50, 5])

p0 = [p0 + 1e-5 * np.random.randn(ndim) for k in range(nwalkers)]

# Set up the sampler.
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob,
    args=(x, y, yerr, bounds))

# Run a burn-in chain and save the final location.
pos, _, _, _ = sampler.run_mcmc(p0, 1000)

# Run the production chain.
sampler.reset()
sampler.run_mcmc(pos, 1000);

# Plots the walkers. Useful for diagnostics, particular of the priors.
# fig, axes = plt.subplots(len(params), 1, sharex=True,
#     figsize=(8, len(params)*3))
# for i,p in enumerate(params):
#     axes[i].plot(sampler.chain[:, :, i].T, color="k", alpha=0.4)
#     axes[i].yaxis.set_major_locator(MaxNLocator(5))
#     axes[i].set_ylabel("${p}$".format(p=p))

# Generates a corner plot (of all variables against each other).
# Useful for diagnostics.
# triangle.corner(sampler.flatchain, bins=50, extents=bounds, labels=params);

samples = sampler.flatchain[:, :2]
m_mcmc, b_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))
print("""FUV - MCMC result:
    m = {m} (+{m_upper}, -{m_lower})
    b = {b} (+{b_upper} -{b_lower})
""".format(m=m_mcmc[0], m_upper=m_mcmc[1], m_lower=m_mcmc[2],
           b=b_mcmc[0], b_upper=b_mcmc[1], b_lower=b_mcmc[2]))

# Setup reference variables for emperical deadtime correction
gcr = np.arange(1,100000,1000)
tec2fdead=5.52e-6
feeclkratio=0.966
dt_old = tec2fdead*(gcr)/feeclkratio
dt_new = 1-(m_mcmc[0] * gcr + b_mcmc[0])/b_mcmc[0]

norm = 0.0
post_prob = np.zeros(len(x))
for i in range(sampler.chain.shape[1]):
    for j in range(sampler.chain.shape[0]):
        ll_fg, ll_bg = sampler.blobs[i][j]
        post_prob += np.exp(ll_fg - np.logaddexp(ll_fg, ll_bg))
        norm += 1
post_prob /= norm

# Plot the prediction.
plt.figure(figsize=(8*1,4*1))
x0 = np.linspace(0,18000, 1000)
A = np.vander(x0, 2)
lines = np.dot(sampler.flatchain[:, :2], A.T)
quantiles = np.percentile(lines, [16, 84], axis=0)
plt.fill_between(x0, quantiles[0], quantiles[1], color="#8d44ad", alpha=0.5)

# Plot the "foreground" points.
ix_bg = np.where(post_prob<.5)
bg_ct = len(ix_bg[0])
plt.errorbar(x[ix_bg], y[ix_bg], yerr=yerr[ix_bg], fmt=",k",
    marker='.',alpha=0.1)

# Plot the "noise" points.
ix_fg = np.where(post_prob>=.5)
fg_ct = len(ix_fg[0])
plt.errorbar(x[ix_fg], y[ix_fg], yerr=yerr[ix_fg], fmt=",k", marker='o', ms=0,
    capsize=0, lw=1, zorder=999, alpha=0.1)
plt.text(5000, 76,
    r'$scr={m}^{{{mp}}}_{{{mm}}}gcr+{b}^{{{bp}}}_{{{bm}}}$'.format(
        m=round(m_mcmc[0],6),b=round(b_mcmc[0],2),
        mp='+{v}'.format(v=round(m_mcmc[1],6)),
        mm='-{v}'.format(v=round(m_mcmc[2],6)),
        bp='+{v}'.format(v=round(b_mcmc[1],2)),
        bm='-{v}'.format(v=round(b_mcmc[2],2))), fontsize=18)

#plt.title('{b} Stim vs. Global Countrate (n={n})'.format(
#    b=band,n=bg_ct+fg_ct),fontsize=16)
plt.xlabel("{b} Global Countrate (ct/s, n={n})".format(
    b=band,n=bg_ct+fg_ct),fontsize=14)
plt.ylabel("{b} Stim Countrate (ct/s)".format(b=band),fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.ylim(68, 78)
plt.xlim(0, 17000)
plt.tight_layout()
plt.savefig('{path}/Fig10a.pdf'.format(path=outpath,b=band),
    format='pdf',dpi=1000)
#plt.close()

# Print a summary
print("""FUV - MCMC result:
    m = {m} (+{m_upper}, -{m_lower})
    b = {b} (+{b_upper} -{b_lower})
""".format(m=m_mcmc[0], m_upper=m_mcmc[1], m_lower=m_mcmc[2],
           b=b_mcmc[0], b_upper=b_mcmc[1], b_lower=b_mcmc[2]))

print 'Tossing ~{p}% of data. ({n} of {m})'.format(
    p=100*bg_ct/(fg_ct+bg_ct),n=bg_ct,m=bg_ct+fg_ct)

#-------------------------------
# NUV linear mixture model
scl = 1.4

inpath = '.'
outpath = '.'

band = 'NUV'
dtdata = {'NUV':pd.read_csv(
    '{path}/NUVStim.csv'.format(path=inpath),index_col=0)}
print dtdata[band].keys()

band = 'NUV'
ix = np.where(np.isfinite(dtdata[band]['globalcount']) &
              np.isfinite(dtdata[band]['stimcount']) &
              (dtdata[band]['exptime']>0))
rawexpt = np.array(dtdata[band]['t1'])[ix]-np.array(dtdata[band]['t0'])[ix]
x = np.array(dtdata[band]['globalcount'])[ix]/rawexpt
y = np.array(dtdata[band]['stimcount'])[ix]/rawexpt
xerr = (np.sqrt(np.array(dtdata[band]['globalcount'])[ix])/
                                        np.array(dtdata[band]['exptime'])[ix])
yerr = (np.sqrt(np.array(dtdata[band]['stimcount'])[ix])/
                                        np.array(dtdata[band]['exptime'])[ix])
print x.min(),x.max(),xerr.min(),xerr.max()
print y.min(),y.max(),yerr.min(),yerr.max()

# Plot the data. Useful for diagnostics.
# plt.figure()
# plt.title('{b} Stim vs. Global Count Rates'.format(b=band))
# plt.xlabel('Global Count Rate (ct/s)')
# plt.ylabel('Stim Count Rate (ct/s)')
# plt.plot(x,y,'.',alpha=0.1,label='{b}'.format(b=band))
#
# plt.figure()
# plt.title('{b} Stim vs. Global Count Rates Errors'.format(b=band))
# plt.xlabel('Global Count Rate Error (ct/s)')
# plt.ylabel('Stim Count Rate Error (ct/s)')
# plt.plot(xerr,yerr,'.',alpha=0.1,label='{b}'.format(b=band))
#
# plt.show()

# Build a linear mixture model. Run MCMC on it.
# Based on http://dan.iel.fm/emcee/current/user/line/
# Define the probabilistic model...
def lnprior(p,bounds):
    # We'll just put reasonable uniform priors on all the parameters.
    if not all([b[0] < v < b[1] for v, b in zip(p, bounds)]):
        return -np.inf
    return 0

# The "foreground" linear likelihood:
def lnlike_fg(p,x,y,yerr):
    m, b, _, M, lnV = p
    model = m * x + b
    return -0.5 * (((model - y) / yerr) ** 2 + 2 * np.log(yerr))

# The "background" outlier likelihood:
def lnlike_bg(p,x,y,yerr):
    _, _, Q, M, lnV = p
    var = np.exp(lnV) + yerr**2
    return -0.5 * ((M - y) ** 2 / var + np.log(var))

# Full probabilistic model.
def lnprob(p,x,y,yerr,bounds):
    m, b, Q, M, lnV = p
    lp = lnprior(p,bounds)
    if not np.isfinite(lp):
        return -np.inf, None
    ll_fg = lnlike_fg(p,x,y,yerr)
    arg1 = ll_fg + np.log(Q)
    ll_bg = lnlike_bg(p,x,y,yerr)
    arg2 = ll_bg + np.log(1.0 - Q)
    ll = np.sum(np.logaddexp(arg1, arg2))
    return lp + ll, (arg1, arg2)

# Initialize the walkers at a reasonable location.
ndim, nwalkers = 5, 32
params = ['m', 'b', 'Q', 'M', 'lnV']
bounds = [(-.001,0), (72, 81), (0, 1), (0, 100), (0, 10)]
p0 = np.array([-0.0005, 79, 0.5, 50, 5])

p0 = [p0 + 1e-5 * np.random.randn(ndim) for k in range(nwalkers)]

# Set up the sampler.
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr, bounds))

# Run a burn-in chain and save the final location.
pos, _, _, _ = sampler.run_mcmc(p0, 500)

# Run the production chain.
sampler.reset()
sampler.run_mcmc(pos, 1000);

# Plot the walkers. Useful for diagnostics.
# fig, axes = plt.subplots(len(params), 1, sharex=True, figsize=(8, len(params)*3))
# for i,p in enumerate(params):
#     axes[i].plot(sampler.chain[:, :, i].T, color="k", alpha=0.4)
#     axes[i].yaxis.set_major_locator(MaxNLocator(5))
#     axes[i].set_ylabel("${p}$".format(p=p))

# Make "triangle" plots of every variable against every other. For diagnostics.
# triangle.corner(sampler.flatchain, bins=50, extents=bounds, labels=params);

samples = sampler.flatchain[:, :2]
m_mcmc, b_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))
print("""MCMC result:
    m = {m} (+{m_upper}, -{m_lower})
    b = {b} (+{b_upper} -{b_lower})
""".format(m=m_mcmc[0], m_upper=m_mcmc[1], m_lower=m_mcmc[2],
           b=b_mcmc[0], b_upper=b_mcmc[1], b_lower=b_mcmc[2]))

gcr = np.arange(1,100000,1000)
tec2fdead=5.52e-6
feeclkratio=0.966

norm = 0.0
post_prob = np.zeros(len(x))
for i in range(sampler.chain.shape[1]):
    for j in range(sampler.chain.shape[0]):
        ll_fg, ll_bg = sampler.blobs[i][j]
        post_prob += np.exp(ll_fg - np.logaddexp(ll_fg, ll_bg))
        norm += 1
post_prob /= norm

plt.figure(figsize=(8*1,4*1))
# Plot the predition.
x0 = np.linspace(0,80000, 1000)
A = np.vander(x0, 2)
lines = np.dot(sampler.flatchain[:, :2], A.T)
quantiles = np.percentile(lines, [16, 84], axis=0)
plt.fill_between(x0, quantiles[0], quantiles[1], color="#8d44ad", alpha=0.5)

# Plot the "bad" points.
ix_bg = np.where(post_prob<.5)
bg_ct = len(ix_bg[0])
plt.errorbar(x[ix_bg], y[ix_bg], yerr=yerr[ix_bg], fmt=",k", marker='.',alpha=0.1)
# Plot the "good" points.
ix_fg = np.where(post_prob>=.5)
fg_ct = len(ix_fg[0])
plt.errorbar(x[ix_fg], y[ix_fg], yerr=yerr[ix_fg], fmt=",k", marker='o', ms=0,
    capsize=0, lw=1, zorder=999, alpha=0.1)
plt.text(15000, 45,
    r'$scr={m}^{{{mp}}}_{{{mm}}}gcr+{b}^{{{bp}}}_{{{bm}}}$'.format(
        m=round(m_mcmc[0],6),b=round(b_mcmc[0],2),
        mp='+{v}'.format(v=round(m_mcmc[1],6)),
        mm='-{v}'.format(v=round(m_mcmc[2],6)),
        bp='+{v}'.format(v=round(b_mcmc[1],2)),
        bm='-{v}'.format(v=round(b_mcmc[2],2))), fontsize=18)

#plt.title('{b} Stim vs. Global Countrate (n={n})'.format(
#    b=band,n=fg_ct+bg_ct),fontsize=16)
plt.xlabel("{b} Global Countrate (ct/s, n={n})".format(
    b=band,n=fg_ct+bg_ct),fontsize=14)
plt.ylabel("{b} Stim Countrate (ct/s)".format(b=band),fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.ylim(40, 75)
plt.xlim(10000, 70000)
plt.tight_layout()
plt.savefig('{path}/Fig10b.pdf'.format(path=outpath,b=band),
    format='pdf',dpi=1000)

# Report some results.
print("""MCMC result:
    m = {m} (+{m_upper}, -{m_lower})
    b = {b} (+{b_upper} -{b_lower})
""".format(m=m_mcmc[0], m_upper=m_mcmc[1], m_lower=m_mcmc[2],
           b=b_mcmc[0], b_upper=b_mcmc[1], b_lower=b_mcmc[2]))

print 'Tossing ~{p}% of data. ({n} of {m})'.format(
    p=100*bg_ct/(fg_ct+bg_ct),n=bg_ct,m=fg_ct+bg_ct)

#################
'''Bin Size Estimates'''

""" Estimated errors (in AB Mag) for a given source brightness (in AB Mag)
and integration depth (in seconds).
"""

import numpy as np
import matplotlib.pyplot as plt
import gPhoton.galextools as gt

outpath = '.'

magrange = np.arange(14,24,1)[::-1]

expt_ratio = 0.8 # Estimate of ration of effective to raw exposure time
t_raw = np.arange(1,1600,1)
t_eff = expt_ratio*t_raw

for sigma in [3]:
    fig = plt.figure(figsize=(8*2,6*2))
    #fig = plt.figure(figsize=(6*2,6*2))
    for i,band in enumerate(['FUV','NUV']):
        plt.subplot(2,1,i+1)
        plt.semilogx()
        plt.xlim(1,1600)
        #plt.title('{b} Bin Depths vs. Error'.format(b=band),fontsize=16)
        plt.xlabel('Exposure Bin Depth (s)',fontsize=14)
        plt.ylabel('{b} {n}{s} Error (AB Mag)'.format(b=band,n=sigma,s=r'$\sigma$'),fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=12)
        for mag in magrange:
            cps = gt.mag2counts(mag,band)
            cps_err = sigma*np.sqrt(cps*t_eff)/t_eff
            mag_err = mag-gt.counts2mag(cps+cps_err,band)
            plt.plot(t_raw,mag_err,label='{m}{u}'.format(m=mag,u=' AB Mag' if mag==max(magrange) else ''))
        for l in [30]:
            plt.axvline(l, color='k', linestyle='dotted', linewidth=2, label='{n} Seconds'.format(n=l))
        plt.legend(fontsize=14)
    plt.savefig('{p}/Fig11.pdf'.format(p=outpath,n=sigma),
        format='pdf',dpi=1000,bbox_inches='tight')
