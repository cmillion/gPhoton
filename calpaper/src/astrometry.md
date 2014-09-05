**Steps:**
1. Use the coadds to find a bunch of bright-ish sources in a deep-ish part of the sky (like CDFS).
2. Compute the astrometric centers for all of those sources from the photons and compare to the coadd centers.
3. Crossmatch the visit level detections to the coadd and compare.

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

%pylab
import gQuery as gq
import galextools as gt
import gphoton_utils as gu
from gAperture import gAperture
band, ra0, dec0, radius = ('FUV', 323.06766667, 0.254, 0.1)
ra0, dec0 = 53.1032558472, -27.7963826072 # PS_CDFS_MOS00
ra0, dec0 = 53.1273118244, -27.8744513656 # CDFS_00
radius = 0.5
maglimit = 20.
aper = 4

out =np.array( gq.getArray(gq.mcat_sources(band,ra0,dec0,0.5,maglimit=maglimit)))
#data = {'ra':out[:,0],'dec':out[:,1],'nmag':out[:,2],'fmag':out[:,3],'fexpt':out[:,9],'nexpt':out[:,10]}
data = {'ra':out[:,0],'dec':out[:,1],'nmag':out[:,2],'fmag':out[:,3]}

skypos = unique_sources(data['ra'],data['dec'],data['fmag'],data['nmag'])

radius = 0.005
plt.ioff()
for pos in skypos:
    print pos
    d = gAperture(band,pos,radius,verbose=2)
    c = get_coadds(band,pos[0],pos[1],radius,maglimit)
    bot,top=gu.model_errors(c['fmag'].mean(),band)
    fig = plt.figure()
    plt.gca().invert_yaxis()
    plt.plot(d['exptime'],d['mag']-gt.apcorrect2(band,radius),'.')
    plt.axis([0,1600.,d['mag'].min()-0.01,d['mag'].max()+0.01])
    plt.plot(top)
    plt.plot(bot)
    plt.savefig(str(pos[0])+'_'+str(pos[1])+'.png')
    for mag in c['fmag']:
        plt.plot(np.arange(1600),np.zeros(1600)+mag-gt.apcorrect1(gt.aper2deg(aper),band))
    plt.close()


