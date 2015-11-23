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
            cols = ['objid','flat_counts','mcat_bg','bg_counts',
                    'flux_bgsub_err','cps_mcatbgsub','counts','mag_mcatbgsub',
                    'cps_err','mag_bgsub','cps_bgsub','detys','flux_bgsub',
                    'flux_err','mag_err_1','cps_bgsub_err','t1_data','bg',
                    'responses','t_mean','cps_mcatbgsub_err','mag_bgsub_err_1',
                    'mag_err_2','t0_data','racent','deccent','mag','exptime',
                    'bg_flat_counts','detxs','t0','t1','mag_mcatbgsub_err_2',
                    'flux','mag_mcatbgsub_err_1','flags','mag_bgsub_err_2',
                    'detrad','cps','flux_mcatbgsub_err','flux_mcatbgsub',
                    'mcat_expt','ra','dec','aper4','aper4_err','mcat_bg']

            spreadsheet = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spreadsheet.writerow(cols)
    return extant_objids

def construct_row(i,band,objid,mcat,data):
    # Note: mcat['skybg'] is in counts per second per square arcseconds
    #       where as gPhoton is reporting cps over the aperture area.
    return (objid,data['flat_counts'][0],data['mcat_bg'][0],
            data['bg_counts'][0],data['flux_bgsub_err'][0],
            data['cps_mcatbgsub'][0],data['counts'][0],data['mag_mcatbgsub'][0],
            data['cps_err'][0],data['mag_bgsub'][0],data['cps_bgsub'][0],
            data['detys'][0],data['flux_bgsub'][0],data['flux_err'][0],
            data['mag_err_1'][0],data['cps_bgsub_err'][0],data['t1_data'][0],
            data['bg'][0],data['responses'][0],data['t_mean'][0],
            data['cps_mcatbgsub_err'][0],data['mag_bgsub_err_1'][0],
            data['mag_err_2'][0],data['t0_data'][0],np.array(data['racent'][0],
            dtype='float32'),np.array(data['deccent'][0],dtype='float32'),
            data['mag'][0],data['exptime'][0],data['bg_flat_counts'][0],
            data['detxs'][0],data['t0'][0],data['t1'][0],
            data['mag_mcatbgsub_err_2'][0],data['flux'][0],
            data['mag_mcatbgsub_err_1'][0],data['flags'][0],
            data['mag_bgsub_err_2'][0],data['detrad'][0],data['cps'][0],
            data['flux_mcatbgsub_err'][0],data['flux_mcatbgsub'][0],
            mcat[band]['expt'][i],np.array(mcat['ra'][i],dtype='float32'),
            np.array(mcat['dec'][i],dtype='float32'),mcat[band][4]['mag'][i],
            mcat[band][4]['err'][i],mcat[band]['skybg'][i])

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
