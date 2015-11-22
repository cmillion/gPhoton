import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import pandas as pd
import gQuery as gq
import galextools as gt
import dbasetools as dt
import MCUtils as mc
from gAperture import gAperture

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
            cols = ['objid','t0','t1','t_raw','t_eff','ra','dec','racent',
                    'deccent','aper4','aper4_err','cps_bgsub','cps',
                    'flux_bgsub','flux',
                    'mag_bgsub','mag','distance','response','skybg',
                    'bg','flags']
            spreadsheet = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spreadsheet.writerow(cols)
    return extant_objids

def construct_row(i,band,objid,mcat,data):
    # Note: mcat['skybg'] is in counts per second per square arcseconds
    #       where as gPhoton is reporting cps over the aperture area.
    return (objid, data['t0'][0], data['t1'][0],
            mcat[band]['expt'][i], data['exptime'][0],
            np.array(mcat['ra'][i],dtype='float32'),
            np.array(mcat['dec'][i],dtype='float32'),
            np.array(data['racent'][0],dtype='float32'),
            np.array(data['deccent'][0],dtype='float32'),
            mcat[band][4]['mag'][i], mcat[band][4]['err'][i],
            data['cps_bgsub'][0], data['cps'][0],
            data['flux_bgsub'][0], data['flux'][0],
            data['mag_bgsub'][0], data['mag'][0],
            mc.distance(data['detxs'],data['detys'],400,400)[0],
            data['responses'][0], mcat[band]['skybg'][i],data['bg'][0],
            data['flags'][0])

def datamaker(band,skypos,outfile,maglimit=20.,margin=0.005,searchradius=0.1,
              radius=gt.aper2deg(4),annulus=[0.0083,0.025],verbose=0):
    extant_objids = file_setup(outfile)
    if extant_objids==False:
        print 'NOT RUNNING!!*!'
        return False
    uniques = dt.find_unique_sources(band,skypos[0],skypos[1],
                                                searchradius,maglimit=maglimit)
    if uniques is None:
        print 'No sources at this position.'
        return
    for pos in uniques:
        mcat = dt.get_mcat_data(pos,margin)
        if not mcat:
            print 'Nothing at {pos}.'.format(pos=pos)
            continue
        extant_objids = file_setup(outfile)
        for i,objid in enumerate(mcat['objid']):
            if objid in extant_objids:
                print 'Already processed.'
                continue
            exp = dt.exp_from_objid(objid)
            if exp[band]['t0']<0:
                print 'skip'
                continue
            print [mcat['ra'][i],mcat['dec'][i]]
            print [exp[band]['t0'],exp[band]['t1']]
            data = gAperture(band,[mcat['ra'][i],mcat['dec'][i]],radius,
                             annulus=annulus,verbose=0,coadd=True,
                             trange=[exp[band]['t0'],exp[band]['t1']],
                             detsize=1.25)
            csv_construct = construct_row(i,band,objid,mcat,data)
            print csv_construct
            with open(outfile,'ab') as csvfile:
                spreadsheet = csv.writer(csvfile, delimiter=',',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
                spreadsheet.writerow(csv_construct)
    return
