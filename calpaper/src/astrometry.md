%pylab
import gQuery as gq
import galextools as gt
import dbasetools as dt
import gphoton_utils as gu
from gAperture import gAperture
band, ra0, dec0, radius = ('FUV', 323.06766667, 0.254, 0.1)
ra0, dec0 = 53.1032558472, -27.7963826072 # PS_CDFS_MOS00
ra0, dec0 = 53.1273118244, -27.8744513656 # CDFS_00
radius = 0.5
maglimit = 20.
aper = 4

data = dt.get_mags(band,ra0,dec0,0.5,maglimit,mode='coadd')
skypos = np.array(dt.parse_unique_sources(data['ra'],data['dec'],
                         data['FUV']['mag'],data['NUV']['mag'],margin=0.001))

# Time range query...
#   select top 10 objid, minPhotoObsDate, maxPhotoObsDate, obs_date, obsdatim, nobs_dat, nobssecs, nobsdati, fobs_dat, fobssecs, fobsdati, nexptime, fexptime
#   from visitphotoobjall as vp
#   inner join imgrun as ir on vp.photoextractid=ir.imgrunid
#   inner join visitphotoextract as vpe on vp.photoextractid=vpe.photoextractid

aper = 4
radius = gt.aper2deg(aper)
annulus = [0.01,0.02]
ac = gt.apcorrect1(radius,band)
plt.ioff()
for i, pos in enumerate(skypos):
    print i, pos
    d = gAperture(band,pos,radius,verbose=2,minexp=30,maxgap=10)#,annulus=annulus)
    c = dt.get_mags(band,pos[0],pos[1],0.0001,maglimit+1,mode='coadd')
    refmag = gt.counts2mag(gt.mag2counts(c[band][aper],band).mean(),band)
    bot,top=gu.model_errors(refmag-ac,band)
    plt.figure()
    plt.gca().invert_yaxis()
    c = dt.get_mags(band,pos[0],pos[1],0.001,maglimit+1,mode='visit')
    plt.plot(c[band]['expt'],c[band][aper]-ac,'x')
    plt.plot(top)
    plt.plot(bot)
    for mag in c[band][aper]:
        plt.plot(np.arange(1600),np.zeros(1600)+mag-ac,color='0.75')
    #plt.plot(d['exptime'],d['mag_bgsub_cheese']-ac,'o')
    #plt.plot(d['exptime'],d['mag_bgsub']-ac,'x')    
    plt.plot(d['exptime'],d['mag']-ac,'.')
    plt.axis([0,1600.,min(d['mag'].min()-ac-0.01,c[band][aper].min()-0.01),
                      max(d['mag'].max()-ac+0.01,c[band][aper].max()+0.01)])
    print "Saving figure."
    plt.savefig(str(i)+'_'+str(pos[0])+'_'+str(pos[1])+'.png')
    plt.close()


