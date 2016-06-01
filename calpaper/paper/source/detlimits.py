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
