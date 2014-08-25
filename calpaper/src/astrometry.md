import gQuery as gq
band, ra0, dec0, radius = ('FUV', 323.06766667, 0.254, 0.1)
ra0, dec0 = 53.1032558472, -27.7963826072 # PS_CDFS_MOS00
ra0, dec0 = 53.1273118244, -27.8744513656 # CDFS_00

out =np.array( gq.getArray(gq.mcat_visit_sources(band,ra0,dec0,0.5)))
data = {'ra':out[:,0],'dec':out[:,1],'nmag':out[:,2],'fmag':out[:,3],'fexpt':out[:,9],'nexpt':out[:,10]}

