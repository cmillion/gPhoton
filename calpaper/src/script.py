%pylab

from regtestutils import *
import galextools as gt
import MCUtils as mc
import pandas as pd
import gphoton_utils as gu
import gFind
from FileUtils import flat_filename

skypos = {'LDS749B':[323.06766667,0.25400000],
          'PS_CDFS_MOS00':[53.1032558472, -27.7963826072],
          'CDFS_00':[53.1273118244, -27.8744513656]}
bands = ['NUV','FUV']

ra0,dec0=skypos['LDS749B']

band = 'FUV'
outfile = 'LDS749B_FUV.csv'
datamaker(band,skypos['LDS749B'],outfile)

band = 'NUV'
outfile = 'LDS749B_NUV.csv'
datamaker(band,skypos['LDS749B'],outfile)

"""To avoid wasting a lot of time on the same tiles and also biasing the
result with the same few super deep tiles (LDS and CDFS, mostly), the
following filters on <5000s coadd depth, more or less, based upon
the most recent coadd level database coverage reference table.
"""
coaddlist = pd.read_csv('coaddList_061714.csv')
for skypos in zip(coaddlist['avaspra'],coaddlist['avaspdec']):
    expt = gFind.gFind(skypos=skypos,band='FUV',quiet=True)['expt'][0]
    if (expt<=5000.) and (expt>0):
        print skypos, expt, True
        datamaker('FUV',skypos,'calrun_FUV-foo.csv',maglimit=24)
    else:
        print skypos, expt, False

coaddlist = pd.read_csv('coaddList_061714.csv')
for skypos in zip(coaddlist['avaspra'],coaddlist['avaspdec']):
    expt = gFind.gFind(skypos=skypos,band='NUV',quiet=True)['expt'][0]
    if (expt<=5000.) and (expt>0):
        print skypos, expt, True
        datamaker('NUV',skypos,'calrun_NUV.csv',maglimit=24)
    else:
        print skypos, expt, False

###############################################################################
def error(data,band,radius,annulus):
    N_a = 1
    N_b0 = (mc.area(annulus[1])-mc.area(annulus[0]))/mc.area(radius)
    N_b = data[band]['bg_eff_area']/mc.area(radius)
    B0 = data[band]['bg']
    B = data[band]['bg_cheese']
    S = gt.mag2counts(data[band]['mag'],band)*data[band]['t_eff']
    s2 = {'bg_cheese_err':(S-B)+(N_a+(N_a**2.)/N_b),
          'bg_err':(S-B0)+(N_a+(N_a**2.)/N_b0)}
    return s2

###############################################################################
# Jake VanderPlas was nice enough to make a clean looking template for us...
from astroML.plotting import setup_text_plots
scl = 1.8
setup_text_plots(fontsize=8*scl, usetex=False)

bands = ['NUV','FUV']
base = 'calrun_'
data = {}
for band in bands:
    data[band] = pd.read_csv('{base}{band}.csv'.format(
                                                        base=base,band=band))
    print '{band} sources: {cnt}'.format(
                                band=band,cnt=data[band]['objid'].shape[0])

"""dMag vs. Mag"""
dmag = {'gphot_cheese':(lambda band:
                        data[band]['aper4']-data[band]['mag_bgsub_cheese']),
        'gphot_nomask':(lambda band:
                        data[band]['aper4']-data[band]['mag_bgsub']),
        'mcat':lambda band: data[band]['aper4']-
            gt.counts2mag(gt.mag2counts(data[band]['mag'],band)-
            data[band]['skybg']*3600**2*mc.area(gt.aper2deg(4)),band)}
for bgmode in dmag.keys():
    for band in bands:
        fig = plt.figure(figsize=(8*scl,4*scl))
        fig.subplots_adjust(left=0.12,right=0.95,wspace=0.02,
                                                        bottom=0.15,top=0.9)
        dmag_err=gu.dmag_errors(100.,band,sigma=1.41)
        # Make a cut on crazy outliers in the MCAT. Also on det radius and expt.
        ix = ((data[band]['aper4']>0) & (data[band]['aper4']<30) &
              (data[band]['distance']<300) & (data[band]['t_eff']<300))
        plt.subplot(1,2,1)
        plt.title('{band} {d}Mag vs. AB Mag (n={n},{bg}_bg)'.format(
                            d=r'$\Delta$',band=band,n=ix.shape[0],bg=bgmode))
        plt.xlabel('AB Magnitude (MCAT)')
        plt.ylabel(r'{d}Magnitude (MCAT-gPhoton)'.format(d=r'$\Delta$'))
        plt.axis([13,23,-1.3,1.3])
        plt.plot(data[band]['aper4'][ix],dmag[bgmode](band)[ix],'.',
                        alpha=0.25,color='k',markersize=5,label='gPhoton Mags')
        plt.plot(dmag_err[0],dmag_err[1],color='k',
                        label='1.41{s}'.format(s=r'$\sigma$'))
        plt.plot(dmag_err[0],dmag_err[2],color='k')
        plt.legend()
        plt.subplot(1,2,2,xticks=[],yticks=[],ylim=[-1.3,1.3])
        plt.title('{d}Magnitude Histogram ({band})'.format(
                                                    d=r'$\Delta$',band=band))
        plt.hist(dmag[bgmode](band)[ix],bins=np.floor(ix.shape[0]/10.),
                        range=[-1.3,1.3],orientation='horizontal',color='k')
        fig.savefig(
            '../calpaper/src/dMag_v_Mag({band},{bg}_bg).png'.format(
                                                        band=band,bg=bgmode))

"""Compare Errors"""
fig = plt.figure(figsize=(4*scl,8*scl))
fig.subplots_adjust(left=0.12,right=0.95,wspace=0.02,bottom=0.05,top=0.95)
for i,band in enumerate(bands):
    s2 = error(data,band,gt.aper2deg(4),[0.0083,0.025])
    dflux = np.sqrt(s2['bg_cheese_err'])/data[band]['t_eff'] # \Delta cps
    flux = gt.mag2counts(data[band]['mag'],band) # cps
    dmag = 2.5*np.log10((flux+dflux)/flux)
    plt.subplot(2,1,i+1,yticks=[],xlim=[0,1])
    plt.title('{band} Magnitude Errors'.format(band=band))
    plt.hist(dmag,bins=100,range=[0,1],color='r',label='gPhoton')
    plt.hist(data[band]['aper4_err'],bins=100,range=[0,1],color='k',
                                                                label='MCAT')
    plt.legend()
fig.savefig('../calpaper/src/mag_err.png')


"""Background Plots
The gPhoton background is scaled to counts in the aperture and the
MCAT skybg is in units of arcsec^2/s, so we'll scale the skybg to the
aperture and put the gPhoton bg in cps.
"""
fig = plt.figure(figsize=(8*scl,4*scl))
fig.subplots_adjust(left=0.12,right=0.95,wspace=0.02,bottom=0.15,top=0.9)
for i,band in enumerate(bands):
    gphot_bg = gt.counts2mag(data[band]['bg_cheese']/data[band]['t_eff'],band)
    mcat_bg = gt.counts2mag(
                    data[band]['skybg']*3600**2*mc.area(gt.aper2deg(4)),band)
    delta = mcat_bg - gphot_bg
    ix = (np.isfinite(delta))
    plt.subplot(1,2,i+1,yticks=[],xlim=[-3,3])
    plt.title('{band} Background {d}Magnitude (MCAT-gPhoton)'.format(band=band,d=r'$\Delta$'))
    plt.hist(delta[ix],bins=500,range=[-3,3],color='k')
fig.savefig('../calpaper/src/dMag_bg.png')

for i,band in enumerate(bands):
    fig = plt.figure(figsize=(8*scl,4*scl))
    fig.subplots_adjust(left=0.12,right=0.95,wspace=0.02,bottom=0.15,top=0.9)
    gphot_bg = gt.counts2mag(data[band]['bg']/data[band]['t_eff'],band)
    mcat_bg = gt.counts2mag(
                    data[band]['skybg']*3600**2*mc.area(gt.aper2deg(4)),band)
    delta = mcat_bg - gphot_bg
    ix = (np.isfinite(delta))
    plt.subplot(1,2,1,xlim=[14,25])
    plt.title('{band} Background vs. Magnitude (MCAT)'.format(band=band))
    plt.plot(data[band]['aper4'][ix],mcat_bg[ix],'.',color='k',alpha=0.2)
    plt.subplot(1,2,2,xlim=[14,25],yticks=[])
    plt.title('{band} Background vs. Magnitude (gPhoton)'.format(band=band))
    plt.plot(data[band]['mag_bgsub_cheese'][ix],gphot_bg[ix],'.',
                                                        color='k',alpha=0.2)
fig.savefig('../calpaper/src/bg_v_mag.png')

"""Astrometry"""
a = 3600
for i,band in enumerate(bands):
    fig = plt.figure(figsize=(8*scl,8*scl))
    fig.subplots_adjust(left=0.12,right=0.95,hspace=0.02,wspace=0.02,
                                                        bottom=0.15,top=0.9)
    dRA = data[band]['ra']-data[band]['racent']
    dDec = data[band]['dec']-data[band]['deccent']
    # dRA v. dDec
    plt.subplot(2,2,1,xticks=[])
    plt.ylabel('{d}RA (arcsec)'.format(d=r'$\Delta$'))
    plt.title('{band} {d}Centroid (MCAT-gPhoton)'.format(
                                                     band=band,d=r'$\Delta$'))
    plt.axis([-0.004*a,0.004*a,-0.004*a,0.004*a])
    plt.plot(dRA*np.cos(data[band]['dec'])*a,dDec*a,'.',alpha=0.1,color='k',
                                                                 markersize=5)
    # dRA
    plt.subplot(2,2,2,yticks=[],xticks=[],ylim=[-0.004*a,0.004*a])
    plt.hist(dRA*np.cos(data[band]['ra'])*a,bins=500,orientation='horizontal',
                                                                    color='k')
    # dDec
    plt.subplot(2,2,3,yticks=[],xlim=[-0.004*a,0.004*a])
    plt.xlabel('{d}Dec (arcsec)'.format(d=r'$\Delta$'))
    plt.gca().invert_yaxis()
    plt.hist(dDec*a,bins=500,color='k')
    fig.savefig('../calpaper/src/dRA_v_dDec({band})'.format(band=band))

fig = plt.figure(figsize=(8*scl,4*scl))
fig.subplots_adjust(left=0.12,right=0.95,wspace=0.02,bottom=0.15,top=0.9)
for i,band in enumerate(bands):
    delta = mc.angularSeparation(data[band]['ra'],data[band]['dec'],
                                 data[band]['racent'],data[band]['deccent'])
    plt.subplot(1,2,i+1,yticks=[],xlim=[0.*a,0.002*a])
    plt.title('{band} Angular Separation (arcsec)'.format(
                                                    band=band,d=r'$\Delta$'))
    plt.hist(delta*a,bins=500,range=[0.*a,0.002*a],color='k')
    fig.savefig('../calpaper/src/angSep({band}).png'.format(band=band))

###############################################################################
"""Deadtime Sanity Checks
According to the calibration paper, the FUV deadtime correction should be
small (~ a few percent), but it is actually bigger than the NUV correction.
"""
fig = plt.figure(figsize=(8,4))
fig.subplots_adjust(left=0.12,right=0.95,wspace=0.02,bottom=0.15,top=0.9)
for i,band in enumerate(bands):
    plt.subplot(1,2,i+1,yticks=[])
    plt.title('{band} Deadtime Ratio Histogram'.format(band=band))
    plt.hist(data[band]['t_eff']/data[band]['t_raw'],bins=100,range=[0.2,1.2])

"""Response Sanity Checks"""
fig = plt.figure(figsize=(8,4))
fig.subplots_adjust(left=0.12,right=0.95,wspace=0.02,bottom=0.15,top=0.9)
for i,band in enumerate(bands):
    flat = mc.get_fits_data(flat_filename(band,'../cal/'))
    plt.subplot(1,2,i+1,yticks=[])
    plt.title('{band} Response Histogram'.format(band=band))
    plt.hist(flat.flatten(),bins=100,range=[0.2,1.2])
    plt.hist(data[band]['response'],bins=100,range=[0.2,1.2])
