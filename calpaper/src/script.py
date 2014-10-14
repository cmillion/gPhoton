from utils import *
import galextools as gt
import MCUtils as mc
import pandas as pd
import gFind

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

plt.figure()
ix = np.where(fuv['aper4']>0)
plt.hist(fuv['aper4'].ix[ix]-fuv['mag_bgsub_cheese'].ix[ix],bins=1000)
ix = np.where(nuv['aper4']>0)
plt.hist(nuv['aper4'].ix[ix]-nuv['mag_bgsub_cheese'].ix[ix],bins=1000)

plt.figure()
plt.title('Delta Dec vs. Delta RA')
plt.plot(fuv['ra'].ix[ix]-fuv['racent'].ix[ix],fuv['dec'].ix[ix]-fuv['deccent'].ix[ix],'.')
plt.plot(nuv['ra'].ix[ix]-nuv['racent'].ix[ix],nuv['dec'].ix[ix]-nuv['deccent'].ix[ix],'x')

plt.figure()
plt.title('SkyBG Delta Mag vs. Mag (within APER4)')
plt.plot(gt.counts2mag(fuv['skybg'].ix[ix]*3600**2*mc.area(gt.aper2deg(4)),'FUV'),gt.counts2mag(fuv['skybg'].ix[ix]*3600**2*mc.area(gt.aper2deg(4)),'FUV')-gt.counts2mag(fuv['bg_cheese'].ix[ix]/fuv['t_eff'].ix[ix],'FUV'),'.')
plt.plot(gt.counts2mag(nuv['skybg'].ix[ix]*3600**2*mc.area(gt.aper2deg(4)),'NUV'),gt.counts2mag(nuv['skybg'].ix[ix]*3600**2*mc.area(gt.aper2deg(4)),'NUV')-gt.counts2mag(nuv['bg_cheese'].ix[ix]/nuv['t_eff'].ix[ix],'NUV'),'x')



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
