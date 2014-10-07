import numpy as np
import matplotlib.pyplot as plt
import csv
import gQuery as gq
import galextools as gt
import dbasetools as dt
import gphoton_utils as gu
import MCUtils as mc
from gAperture import gAperture

def get_mcat_data(pos,rad):
    out = np.array(gq.getArray(gq.mcat_visit_sources(pos[0],pos[1],rad)))
    return {'objid':np.array(out[:,0],dtype='int64'),
            'ra':np.array(out[:,1],dtype='float32'),
            'dec':np.array(out[:,2],dtype='float32'),
            'NUV':{'mag':np.array(out[:,3],dtype='float32'),
                   'skybg':np.array(out[:,5],dtype='float32'),
                   'expt':np.array(out[:,9],dtype='float32'),
                   1:np.array(out[:,19],dtype='float32'),
                   2:np.array(out[:,20],dtype='float32'),
                   3:np.array(out[:,21],dtype='float32'),
                   4:np.array(out[:,22],dtype='float32'),
                   5:np.array(out[:,23],dtype='float32'),
                   6:np.array(out[:,24],dtype='float32'),
                   7:np.array(out[:,25],dtype='float32')},
            'FUV':{'mag':np.array(out[:,4],dtype='float32'),
                   'skybg':np.array(out[:,6],dtype='float32'),
                   'expt':np.array(out[:,10],dtype='float32'),
                   1:np.array(out[:,12],dtype='float32'),
                   2:np.array(out[:,13],dtype='float32'),
                   3:np.array(out[:,14],dtype='float32'),
                   4:np.array(out[:,15],dtype='float32'),
                   5:np.array(out[:,16],dtype='float32'),
                   6:np.array(out[:,17],dtype='float32'),
                   7:np.array(out[:,18],dtype='float32') } }

def file_setup(outfile):
    extant_objids = []
    if os.path.exists(outfile):
        print 'This file exists.'
        try:
            extant_objids=np.array(pd.read_csv(outfile)['objid']).tolist()
        except:
            print 'And nonstandard!'
            # Raise an exception?
            return False
    else:
        # Initialize the file with a header
        with open(outfile, 'wb') as csvfile:
            cols = ['objid','t0','t1','t_raw','t_eff','ra','dec','racent','deccent','aper4','mag_bgsub_cheese','mag_bgsub','mag','distance','response','skybg','bg','bg_cheese']
            spreadsheet = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spreadsheet.writerow(cols)
    return extant_objids

def construct_row(band,objid,mcat,data):
    # Note: mcat['skybg'] is in counts per second per square arcseconds
    #       where as gPhoton is reporting cps over the aperture area.
    return (objid, data['t0'][0], data['t1'][0], 
            mcat[band]['expt'][i], data['exptime'][0],
            mcat['ra'][i], mcat['dec'][i],
            data['racent'][0], data['deccent'][0],
            mcat[band][4][i]+18.82, data['mag_bgsub_cheese'][0],
            data['mag_bgsub'][0], data['mag'][0],
            mc.distance(data['detxs'],data['detys'],400,400)[0],
            data['responses'][0], mcat[band]['skybg'][i],
            data['bg']['simple'][0], data['bg']['cheese'][0])

def datamaker(band,skypos,outfile,maglimit=20.,detsize=0.5,
                                radius=gt.aper2deg(4),annulus=[0.01,0.02]):
    extant_objids = file_setup(outfile)
    if extant_objids==False:
        print 'NOT RUNNING!!*!'
        return False
    uniques = dt.find_unique_sources(band,skypos[0],skypos[1],
                                                    detsize,maglimit=maglimit)
    for pos in uniques:
        mcat = get_mcat_data(pos,0.005)
        for i,objid in enumerate(mcat['objid']):
            if objid in extant_objids:
                print 'Already processed.'
                continue
            exp = dt.exp_from_objid(objid)
            if exp[band]['t0']<0:
                print 'skip'
                continue
            data = gAperture(band,[mcat['ra'][i],mcat['dec'][i]],radius,annulus=annulus,verbose=0,trange=[exp[band]['t0'],exp[band]['t1']],coadd=True)
            if data['mag_bgsub_cheese'] and np.isfinite(data['mag_bgsub_cheese']):
                csv_construct = construct_row(band,objid,mcat,data)
                print csv_construct
                with open(outfile,'ab') as csvfile:
                    spreadsheet = csv.writer(csvfile, delimiter=',',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    spreadsheet.writerow(csv_construct)
            else:
                print 'no exp'
    return

skypos = {'LDS749B':[323.06766667,0.25400000],
          'PS_CDFS_MOS00':[53.1032558472, -27.7963826072],
          'CDFS_00':[53.1273118244, -27.8744513656]}
bands = ['NUV','FUV']

ra0,dec0=skypos['LDS749B']
band = bands[1]
outfile = 'test.csv'
datamaker(bands[1],skypos['LDS749B'],outfile)


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
