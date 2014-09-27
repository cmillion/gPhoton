import numpy as np
import matplotlib.pyplot as plt
import gQuery as gq
import galextools as gt
import dbasetools as dt
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

def exp_from_objid(objid):
    out = np.array(gq.getArray(gq.mcat_objid_search(objid)))
    return {'NUV':{'expt':np.array(out[:,7],dtype='float')[0],'t0':np.array(out[:,9],dtype='float64')[0]-gt.gpssecs,'t1':np.array(out[:,10],dtype='float64')[0]-gt.gpssecs},'FUV':{'expt':np.array(out[:,8],dtype='float')[0],'t0':np.array(out[:,11],dtype='float64')[0]-gt.gpssecs,'t1':np.array(out[:,12],dtype='float64')[0]-gt.gpssecs}}

skypos = {'LDS749B':[323.06766667,0.25400000],
          'PS_CDFS_MOS00':[53.1032558472, -27.7963826072],
          'CDFS_00':[53.1273118244, -27.8744513656]}
band = ['NUV','FUV']

maglimit = 23.
coadds = dt.get_mags(band,skypos['CDFS_00'][0],skypos['CDFS_00'][1],0.5,
                                                        maglimit,mode='coadd')
uniques = np.array(dt.parse_unique_sources(coadds['ra'],coadds['dec'],
                       coadds['FUV']['mag'],coadds['NUV']['mag'],margin=0.001))

# LDS749B
src = 'LDS749B'
refmag = {'NUV':14.71,'FUV':15.57}
radius = gt.aper2deg(4)
annulus = [0.01,0.02]

#data = gAperture(band[1],skypos[src],0.005,annulus=annulus,verbose=2,
#                 minexp=30,maxgap=10)
#data['dist']=np.sqrt((data['detxs']-400.)**2. + (data['detys']-400.)**2.)

# Grab the per-visit MCAT magnitudes

out = np.array(gq.getArray(gq.mcat_visit_sources(skypos[src][0],
                                                 skypos[src][1],0.005)))
mcat = {'objid':out[:,0],'ra':out[:,1],'dec':out[:,2],'NUV':{'mag':out[:,3],'expt':out[:,9],1:out[:,19],2:out[:,20],3:out[:,21],4:out[:,22],5:out[:,23],6:out[:,24],7:out[:25]},'FUV':{'mag':out[:,4],'expt':out[:,10],1:out[:,12],2:out[:,13],3:out[:,14],4:out[:,15],5:out[:,16],6:out[:,17],7:out[:,18]}}

mag = np.array([])
dmag = np.array([])
for i,objid in enumerate(mcat['objid']):
    out = exp_from_objid(objid)
    if out[band[1]]['t0']<0:
        continue
    data = gAperture(band[1],skypos[src],radius,annulus=annulus,verbose=0,trange=[out[band[1]]['t0'],out[band[1]]['t1']],coadd=True)
    if data['mag_bgsub_cheese']:
        print objid,[out[band[1]]['t0'],out[band[1]]['t1']],data['mag_bgsub_cheese'],float(mcat[band[1]][4][i])+18.82,float(mcat[band[1]][4][i])+18.82-data['mag_bgsub_cheese']
        dmag = np.append(dmag,(float(mcat[band[1]][4][i])+18.82-data['mag_bgsub_cheese'])[0])
        mag = np.append(mag,float(mcat[band[1]][4][i])+18.82)

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






