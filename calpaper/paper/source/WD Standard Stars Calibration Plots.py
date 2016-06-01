
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
legs = {}
legs['FUV'] = np.array([gq.getArray(
            gq.obstype(o))[0][5] for o in np.array(data['FUV']['objid'])])

# Provide some baseline statistics for the sample
ix = {}
for band in bands:
    # The bitwise_and ignores the annulus mask flags because using MCAT bg.
    ix[band] = np.where((data[band][apertxt]>0) &
                        (data[band]['detrad']<200) &
                        (data[band]['t0']<961986575.) &
        (np.bitwise_and(np.array(data[band]['flags'].values,dtype='int16'),
                                        0b00111111)==0))# & (legs[band]>3)
    print '{b}: {m} / {n}'.format(b=band,n=len(data[band][apertxt]),
        m=len(ix[band][0]))
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
plt.plot(np.array(legs[band])[ix[band]],
                        np.array(dmag)[ix[band]],'.',alpha=0.3,color='k')
plt.subplot(1,2,2,yticks=[])
plt.hist(np.array(dmag)[ix[band]],bins=50,range=dmagrange,
                 orientation='horizontal',color='k',histtype='step',normed=1)
plt.axhline(-0.015, color='g', linestyle='solid', linewidth=1)
plt.axhline(-0.05, color='g', linestyle='solid', linewidth=1)
plt.savefig('{path}/Fig09a.pdf'.format(path=outpath),
    format='pdf',dpi=1000,bbox_inches='tight')

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
                format='pdf',dpi=1000,bbox_inches='tight')

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
    plt.savefig('{path}/Fig08a.pdf'.format(path=outpath),
        format='pdf',dpi=1000,bbox_inches='tight')

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
                format='pdf',dpi=1000,bbox_inches='tight')

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
    plt.savefig('{path}/Fig08b.pdf'.format(path=outpath),
        format='pdf',dpi=1000,bbox_inches='tight')

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
        format='pdf',dpi=1000,bbox_inches='tight')

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
        format='pdf',dpi=1000,bbox_inches='tight')

#plt.close('all')
