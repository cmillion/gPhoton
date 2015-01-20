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
                    'deccent','aper4','aper4_err','mag_bgsub_cheese',
                    'mag_bgsub','mag','distance','response','skybg',
                    'bg','bg_cheese','bg_eff_area', 'bg_sigmaclip',
                    'mag_bgsub_sigmaclip']
            spreadsheet = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spreadsheet.writerow(cols)
    return extant_objids

def construct_row(i,band,objid,mcat,data):
    # Note: mcat['skybg'] is in counts per second per square arcseconds
    #       where as gPhoton is reporting cps over the aperture area.
    return (objid, data['t0'][0], data['t1'][0],
            mcat[band]['expt'][i], data['exptime'][0],
            mcat['ra'][i], mcat['dec'][i],
            data['racent'][0], data['deccent'][0],
            mcat[band][4]['mag'][i], mcat[band][4]['err'][i],
            data['mag_bgsub_cheese'][0],
            data['mag_bgsub'][0], data['mag'][0],
            mc.distance(data['detxs'],data['detys'],400,400)[0],
            data['responses'][0], mcat[band]['skybg'][i],
            data['bg']['simple'][0], data['bg']['cheese'][0],
            data['bg']['eff_area'], data['bg']['sigmaclip'][0],
            data['mag_bgsub_sigmaclip'][0])

def datamaker(band,skypos,outfile,maglimit=20.,detsize=0.5,
                                radius=gt.aper2deg(4),annulus=[0.0083,0.025]):
    """Note: If you wanted to change the default annulus, then a good starting
    point would be [0.0083,0.025] (i.e. 30" to 90").
    """
    extant_objids = file_setup(outfile)
    if extant_objids==False:
        print 'NOT RUNNING!!*!'
        return False
    uniques = dt.find_unique_sources(band,skypos[0],skypos[1],
                                                    detsize,maglimit=maglimit)
    for pos in uniques:
        mcat = dt.get_mcat_data(pos,0.005)
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
            data = gAperture(band,[mcat['ra'][i],mcat['dec'][i]],radius,
                             annulus=annulus,verbose=0,coadd=True,
                             trange=[exp[band]['t0'],exp[band]['t1']])
            if (data['mag_bgsub_cheese'] and
                                        np.isfinite(data['mag_bgsub_cheese'])):
                csv_construct = construct_row(i,band,objid,mcat,data)
                print csv_construct
                with open(outfile,'ab') as csvfile:
                    spreadsheet = csv.writer(csvfile, delimiter=',',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    spreadsheet.writerow(csv_construct)
            else:
                print 'no exp'
    return
