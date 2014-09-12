import numpy as np
import matplotlib.pyplot as plt
import gQuery as gq
import galextools as gt
import gphoton_utils as gu
from gAperture import gAperture

def distance(a,b,c,d):
    return np.sqrt( (a-c)**2. + (b-d)**2. )

def unique_sources(ras,decs,fmags,nmags,margin=0.005):
    skypos = zip(ras,decs)
    for i,pos in enumerate(skypos):
        ix = np.where(distance(pos[0],pos[1],ras,decs)<=margin)
        skypos[i]=[ras[ix].mean(),decs[ix].mean()]
    a = skypos #unique_sources(data['ra'],data['dec'])
    b = []
    for i in a:
        if not (i in b):
            b+=[i]
    return b

def get_coadds(band,ra0,dec0,radius,maglimit):
    out =np.array( gq.getArray(gq.mcat_sources(band,ra0,dec0,radius,maglimit=maglimit)))
    return {'ra':out[:,0],'dec':out[:,1],'nmag':out[:,2],'fmag':out[:,3]}

skypos = {'LDS749B':[323.06766667,0.25400000],
          'PS_CDFS_MOS00':[53.1032558472, -27.7963826072],
          'CDFS_00':[53.1273118244, -27.8744513656]}
band = ['NUV','FUV']

# LDS749B
src = 'LDS749B'
refmag = {'NUV':14.71,'FUV':15.57}
radius = gt.aper2deg(4)
annulus = [0.01,0.02]

data = gAperture(band[1],skypos[src],0.005,annulus=annulus,verbose=2,
                 minexp=30,maxgap=10)

# Grab the per-visit MCAT magnitudes
out = np.array(gq.getArray(gq.mcat_visit_sources(band[1],skypos[src][0],
                                                 skypos[src][1],0.005)))
mcat = {'ra':out[:,0],'dec':out[:,1],'NUV':{'mag':out[:,2],'expt':out[:,8],1:out[:,18],2:out[:,19],3:out[:,20],4:out[:,21],5:out[:,22],6:out[:,23],7:out[:24]},'FUV':{'mag':out[:,3],'expt':out[:,9],1:out[:,11],2:out[:,12],3:out[:,13],4:out[:,14],5:out[:,15],6:out[:,16],7:out[:,17]}}

data['dist']=np.sqrt((data['detxs']-400.)**2. + (data['detys']-400.)**2.)

# Set up the plot
plt.figure()
plt.gca().invert_yaxis()
plt.title(str(src)+' - '+str(band[1]))
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
plt.plot(mcat['fexpt'][ix],mcat['fmag'][ix]-gt.apcorrect1(radius,band[1]),'.')
# Plot the gAperture data
ix = np.where(data['dist']<200)
plt.plot(data['exptime'][ix],data['mag'][ix]-gt.apcorrect1(0.005,band[1]),'o')
plt.plot(data['exptime'],data['mag']-gt.apcorrect1(0.005,band[1]),'x')






