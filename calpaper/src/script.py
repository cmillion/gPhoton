from utils import *
import galextools as gt
import MCUtils as mc
import pandas as pd
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

# Let's do this on more reasonable data...
coaddlist = pd.read_csv('coaddList_061714.csv')
for skypos in zip(coaddlist['avaspra'],coaddlist['avaspdec']):
    expt = gFind.gFind(skypos=skypos,band='FUV',quiet=True)['expt'][0]
    if (expt<=5000.) and (expt>0):
        print skypos, expt, True
        datamaker('FUV',skypos,'calrun_FUV.csv')
    else:
        print skypos, expt, False


coaddlist = pd.read_csv('coaddList_061714.csv')
for skypos in zip(coaddlist['avaspra'],coaddlist['avaspdec']):
    expt = gFind.gFind(skypos=skypos,band='NUV',quiet=True)['expt'][0]
    if (expt<=5000.) and (expt>0):
        print skypos, expt, True
        datamaker('NUV',skypos,'calrun_NUV.csv')
    else:
        print skypos, expt, False


fuv=pd.read_csv('calrun_FUV.csv')
nuv=pd.read_csv('calrun_NUV.csv')
print 'FUV samples: {cnt} (blue)'.format(cnt=fuv.shape[0])
print 'NUV samples: {cnt} (red)'.format(cnt=nuv.shape[0])

# NOTE: There is a ~15% systematic offset in FUV!
plt.figure()
plt.title('FUV Delta Mag vs. Mag')
dmag = fuv['aper4']-fuv['mag_bgsub_cheese']
ix = np.where((fuv['aper4']>0) & (fuv['aper4']<30))
plt.axis([10,24,-1.3,1.3])
plt.plot(fuv['aper4'].ix[ix],dmag.ix[ix],'.',alpha=0.1)
plt.figure()
plt.title('FUV Delta Mag Histogram')
plt.axis([-1.3,1.3,0,425])
plt.hist(dmag.ix[ix],bins=500,range=[-1.3,1.3])

plt.figure()
plt.title('NUV Delta Mag vs. Mag')
dmag = nuv['aper4']-nuv['mag_bgsub_cheese']
ix = np.where((nuv['aper4']>0) & (nuv['aper4']<30))
plt.plot(nuv['aper4'].ix[ix],dmag.ix[ix],'.',alpha=0.1)
plt.figure()
plt.title('NUV Delta Mag Histogram')
plt.axis([-1.3,1.3,0,4500])
plt.hist(dmag.ix[ix],bins=500,range=[-1.3,1.3])

# According to the calibration paper, the FUV deadtime correction should be
# small (~ a few percent), but it is actually bigger than the NUV correction.
plt.figure()
plt.title('FUV Deadtime Ratio Histogram')
plt.hist(fuv['t_eff']/fuv['t_raw'],bins=100,range=[0.2,1.2])
plt.figure()
plt.title('NUV Deadtime Ratio Histogram')
plt.hist(nuv['t_eff']/nuv['t_raw'],bins=100,range=[0.2,1.2])

# So what happens to Delta Mag if we estimate ~3% FUV deadtime?
dtime = 0.03
plt.figure()
plt.title('FUV Delta Mag Histogram ({dt}% deadtime)'.format(dt=dtime*100))
plt.hist(fuv['aper4']-gt.counts2mag(
  gt.mag2counts(fuv['mag'],'FUV')*fuv['t_eff']/(fuv['t_raw']*(1-dtime)),'FUV'),
                                                        range=[-1,1],bins=100)

# Do a sanity check of the reponse sampling
flat_FUV = mc.get_fits_data(flat_filename('FUV','../cal/'))
flat_NUV = mc.get_fits_data(flat_filename('NUV','../cal/'))
plt.figure()
plt.hist(flat_FUV.flatten(),bins=100,range=[0.2,1.2])
plt.hist(fuv['response'],bins=100,range=[0.2,1.2])
plt.figure()
plt.hist(flat_NUV.flatten(),bins=100,range=[0.2,1.2])
plt.hist(nuv['response'],bins=100,range=[0.2,1.2])


# RA v Dec
plt.figure()
plt.title('Delta Dec vs. Delta RA')
plt.plot(fuv['ra'].ix[ix]-fuv['racent'].ix[ix],fuv['dec'].ix[ix]-fuv['deccent'].ix[ix],'.')
plt.plot(nuv['ra'].ix[ix]-nuv['racent'].ix[ix],nuv['dec'].ix[ix]-nuv['deccent'].ix[ix],'x')

# Background plots
gphot_back = gt.counts2mag(fuv['bg_cheese']/fuv['t_eff'],'FUV')
mcat_back = gt.counts2mag(fuv['skybg']*3600**2*mc.area(gt.aper2deg(4)),'FUV')
delta = mcat_back - gphot_back
plt.figure()
plt.title('gPhoton Background Histogram')
ix = np.where(np.isfinite(gphot_back))
plt.axis([18,26,0,300])
plt.hist(gphot_back.ix[ix],bins=500,range=[18,26])
plt.figure()
plt.title('MCAT Background Histogram')
ix = np.where(np.isfinite(mcat_back))
plt.axis([18,26,0,300])
plt.hist(mcat_back.ix[ix],bins=500,range=[18,26])
plt.figure()
plt.title('MCAT-gPhoton Background Histogram')
ix = np.where(np.isfinite(delta))
plt.axis([-3,3,0,350])
plt.hist(delta.ix[ix],bins=500,range=[-3,3])





plt.plot(gt.counts2mag(fuv['skybg']*3600**2*mc.area(gt.aper2deg(4)),'FUV'),delta,'.',alpha=0.1)


plt.plot(gt.counts2mag(nuv['skybg'].ix[ix]*3600**2*mc.area(gt.aper2deg(4)),'NUV'),
         gt.counts2mag(nuv['skybg'].ix[ix]*3600**2*mc.area(gt.aper2deg(4)),'NUV')
        -gt.counts2mag(nuv['bg_cheese'].ix[ix]/nuv['t_eff'].ix[ix],'NUV'),'x')



#######################################


# Set up the plot
plt.figure()
plt.gca().invert_yaxis()
plt.title(str(src)+' - '+str(band))
plt.xlabel('Exposure time (s)')
plt.ylabel('AB Magnitude')
# Make the y dimensions a little bigger than the data
plt.axis([0,1600.,data['mag'].min()-0.01,data['mag'].max()+0.01])
# Plot the reference magnitude along with upper and lower bounds vs. expt
top,bot=gu.model_errors(refmag['FUV'],band,sigma=5)
plt.plot(top)
plt.plot(bot)
plt.plot(np.ones(1600)*refmag['FUV'])
# Plot the MCAT data
ix = np.where(mcat['fmag']>0)


#####################################
# Set up the plot
plt.figure()
plt.gca().invert_yaxis()
plt.title(str(src)+' - '+str(band))
plt.xlabel('Exposure time (s)')
plt.ylabel('AB Magnitude')
# Make the y dimensions a little bigger than the data
plt.axis([0,1600.,data['mag'].min()-0.01,data['mag'].max()+0.01])
