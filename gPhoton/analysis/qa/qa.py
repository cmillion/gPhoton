'''
Plot:
(1) cps vs. t w/ 1/3/5 sigma errors
(2) responses vs. t
(3) detrad vs. t
(4) detx vs. t
(5) dety vs. t
(6..) parse out the flags
'''

from gPhoton.gphoton_utils import read_lc
import gPhoton.dbasetools as dt
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import gPhoton.gphoton_utils as gu
import os

def crosscorr_title(a,b):
    return '{pearsonr}, {spearmanr}, {kendalltau}'.format(
        pearsonr=round(scipy.stats.pearsonr(a,b)[0],2),
        spearmanr=round(scipy.stats.spearmanr(a,b)[0],2),
        kendalltau=round(scipy.stats.kendalltau(a,b)[0],2))

def checkplot(csvfile,outfile=None,format='png',maxgap=500,imscale=4,
              nplots=10,cleanup=False):
    lc = read_lc(csvfile)
    tranges = dt.distinct_tranges(np.array(lc['t0']),maxgap=500)
    stepsz = np.median(lc['t1']-lc['t0']) # sort of a guess at stepsz
    n=2 # temporary hacking variable for number of rows in figure
    for j in range(np.int(np.ceil(len(tranges)/float(nplots)))):
        plt.figure(
            figsize=(imscale*(len(tranges[j*nplots:(j+1)*nplots])),imscale*n))
        for i, trange in enumerate(tranges[j*nplots:(j+1)*nplots]):
            # Countrate
            plt.subplot(n, len(tranges[j*nplots:(j+1)*nplots]),i+1+len(tranges[j*nplots:(j+1)*nplots])*0)
            plt.xticks([])
            if i==0:
                plt.ylabel('cps_mcatbgsub')
            plt.ylim(
                np.array(lc['cps_mcatbgsub']-2*5*lc['cps_mcatbgsub_err']).min(),
                np.array(lc['cps_mcatbgsub']+2*5*lc['cps_mcatbgsub_err']).max())
            time_ix = np.where((np.array(lc['t0'])>=trange[0]) &
                                           (np.array(lc['t1'])<=trange[1]))
            if not len(time_ix[0]):
                continue
            tlim = (np.array(lc['t0'])[time_ix].min()-stepsz,
                    np.array(lc['t1'])[time_ix].max()+stepsz)
            plt.xlim(tlim[0],tlim[1])
            for nsigma in [5]:
                plt.errorbar(np.array(lc['t_mean'])[time_ix],
                    np.array(lc['cps_mcatbgsub'])[time_ix],
                    yerr=nsigma*np.array(lc['cps_mcatbgsub_err'])[time_ix],fmt='k.')
            flag_ix = np.where(np.array(lc['flags'])[time_ix]>0)
            plt.plot(np.array(lc['t_mean'])[time_ix][flag_ix],
                np.array(lc['cps_mcatbgsub'])[time_ix][flag_ix],'rx')

            # Detector Radius
            plt.subplot(n, len(tranges[j*nplots:(j+1)*nplots]),i+1+len(tranges[j*nplots:(j+1)*nplots])*1)
            plt.xticks([])
            if i==0:
                plt.ylabel('detrad')
            plt.xlim(tlim[0],tlim[1])
            plt.ylim(np.array(lc['detrad'])[time_ix].min()-2,
                     np.array(lc['detrad'])[time_ix].max()+2)
            plt.plot(np.array(lc['t_mean'])[time_ix],
                     np.array(lc['detrad'])[time_ix],'k.')
            plt.plot(np.array(lc['t_mean'])[time_ix][flag_ix],
                     np.array(lc['detrad'])[time_ix][flag_ix],'rx')
            plt.title(crosscorr_title(
                        np.array(lc['cps_mcatbgsub_err'])[time_ix],
                        np.array(lc['detrad'])[time_ix]))

            # # Response
            # plt.subplot(n, len(tranges),i+1+len(tranges)*2)
            # plt.xticks([])
            # if i==0:
            #     plt.ylabel('responses')
            # plt.xlim(tlim[0],tlim[1])
            # plt.ylim(np.array(lc['responses'])[time_ix].min()-.1,
            #          np.array(lc['responses'])[time_ix].max()+.1)
            # plt.plot(np.array(lc['t_mean'])[time_ix],
            #          np.array(lc['responses'])[time_ix],'k.')
            # plt.title(crosscorr_title(
            #             np.array(lc['cps_mcatbgsub_err'])[time_ix],
            #             np.array(lc['responses'])[time_ix]))
            # plt.plot(np.array(lc['t_mean'])[time_ix][flag_ix],
            #     np.array(lc['responses'])[time_ix][flag_ix],'rx')
            #
            # # Exposure Time
            # plt.subplot(n, len(tranges),i+1+len(tranges)*3)
            # plt.xticks([])
            # if i==0:
            #     plt.ylabel('exptime')
            # plt.xlim(tlim[0],tlim[1])
            # plt.ylim(np.array(lc['exptime'])[time_ix].min()-stepsz/2.,
            #          np.array(lc['exptime'])[time_ix].max()+stepsz/2.)
            # plt.plot(np.array(lc['t_mean'])[time_ix],
            #          np.array(lc['exptime'])[time_ix],'k.')
            # plt.plot(np.array(lc['t_mean'])[time_ix][flag_ix],
            #     np.array(lc['exptime'])[time_ix][flag_ix],'rx')
            # plt.title(crosscorr_title(
            #             np.array(lc['cps_mcatbgsub_err'])[time_ix],
            #             np.array(lc['exptime'])[time_ix]))

        plt.tight_layout()
        if outfile:
            if len(tranges)>nplots:
                plt.savefig('{base}_{j}{type}'.format(
                    base=outfile[:-4],j=j,type=outfile[-4:]),dpi=300)
            else:
                plt.savefig(outfile,dpi=300)
        if cleanup:
            plt.close('all')
    return

def isvariable(csvfile,nsigma=5,spearmanr_cutoff=0.45,nbin_cutoff=5):
    variable, n = False, 0
    if not os.path.exists(csvfile):
        print 'No data: {f}'.format(f=csvfile)
        return False
    out = gu.read_lc(csvfile)
    tranges = dt.distinct_tranges(np.sort(out['t0']),maxgap=500)
    for trange in tranges:
        time_ix = np.where((np.array(out['t0'])>=trange[0]) &
                                            (np.array(out['t1'])<=trange[1]))
        flag_ix = np.where(np.array(out['flags'])[time_ix]==0)
        if len(flag_ix[0])<nbin_cutoff:
            continue # Not enough bins to do anything interesting
        corrcoeff = scipy.stats.spearmanr(
                                    np.array(out['cps_mcatbgsub'])[time_ix],
                                    np.array(out['detrad'])[time_ix])[0]
        if np.abs(corrcoeff)>spearmanr_cutoff:
            continue # Strong correllation to detector radius

        upper = np.array(out['cps_mcatbgsub']+
                            nsigma*out['cps_mcatbgsub_err'])[time_ix][flag_ix]
        lower = np.array(out['cps_mcatbgsub']-
                            nsigma*out['cps_mcatbgsub_err'])[time_ix][flag_ix]
        for u in upper:
            for l in lower:
                if l>u:
                    n+=1
                    if n>1: # Must be more than one significant bin.
                        variable = True
    return variable
