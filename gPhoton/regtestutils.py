"""
.. module:: regtestutils
   :synopsis: Functions for end-to-end photometric regression testing with
       emphasis on comparing gAperture values to MCAT values.

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

from __future__ import absolute_import, division, print_function
# Core and Third Party imports.
import csv
import os
import numpy as np
import pandas as pd
# gPhoton imports.
import gPhoton.dbasetools as dt
import gPhoton.galextools as gt
from gPhoton import gAperture

# ------------------------------------------------------------------------------
def file_setup(outfile):
    """
    Checks for a CSV file in which to put all the data and initializes it
        if it hasn't already been created. Loads already processed data, if any,
        in order to continue interrupted runs.

    :param outfile: Name of output file to make.

    :type outfile: str

    :returns: numpy.ndarray -- Set of object IDs that are already processed.
    """

    extant_objids = []

    if os.path.exists(outfile):
        print('This file exists.')
        try:
            extant_objids = np.array(pd.read_csv(outfile)['objid']).tolist()
        except:
            print('And nonstandard!')
            # Raise an exception?
            return False
    else:
        # Initialize the file with a header
        with open(outfile, 'wb') as csvfile:
            cols = ['objid', 'flat_counts', 'mcat_bg', 'bg_counts',
                    'flux_bgsub_err', 'cps_mcatbgsub', 'counts',
                    'mag_mcatbgsub', 'cps_err', 'mag_bgsub', 'cps_bgsub',
                    'detys', 'flux_bgsub', 'flux_err', 'mag_err_1',
                    'cps_bgsub_err', 't1_data', 'bg', 'responses', 't_mean',
                    'cps_mcatbgsub_err', 'mag_bgsub_err_1', 'mag_err_2',
                    't0_data', 'racent', 'deccent', 'mag', 'exptime',
                    'bg_flat_counts', 'detxs', 't0', 't1',
                    'mag_mcatbgsub_err_2', 'flux', 'mag_mcatbgsub_err_1',
                    'flags', 'mag_bgsub_err_2', 'detrad', 'cps',
                    'flux_mcatbgsub_err', 'flux_mcatbgsub', 'mcat_expt', 'ra',
                    'dec', 'aper4', 'aper4_err', 'mcat_bg',
                    'aper7', 'aper7_err']

            spreadsheet = csv.writer(csvfile, delimiter=',', quotechar='|',
                                     quoting=csv.QUOTE_MINIMAL)
            spreadsheet.writerow(cols)

    return extant_objids
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def construct_row(i, band, objid, mcat, data):
    """
    Assemble gAperture and MCAT data into a CSV row.

    :param i: The index of the row to collect values for.

    :type i: int

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param objid: The GALEX MCAT object ID.

    :type objid: long

    :param mcat: Object containing MCAT data.

    :type mcat: dict

    :param data: Object containing gAperture data.

    :type data: dict

    :returns: tuple -- The CSV row to output.
    """

    # Note: mcat['skybg'] is in counts per second per square arcseconds
    # where as gPhoton is reporting cps over the aperture area.
    return (objid, data['flat_counts'][0], data['mcat_bg'][0],
            data['bg_counts'][0], data['flux_bgsub_err'][0],
            data['cps_mcatbgsub'][0], data['counts'][0],
            data['mag_mcatbgsub'][0], data['cps_err'][0], data['mag_bgsub'][0],
            data['cps_bgsub'][0], data['detys'][0], data['flux_bgsub'][0],
            data['flux_err'][0], data['mag_err_1'][0], data['cps_bgsub_err'][0],
            data['t1_data'][0], data['bg'][0], data['responses'][0],
            data['t_mean'][0], data['cps_mcatbgsub_err'][0],
            data['mag_bgsub_err_1'][0], data['mag_err_2'][0],
            data['t0_data'][0], np.array(data['racent'][0], dtype='float32'),
            np.array(data['deccent'][0], dtype='float32'), data['mag'][0],
            data['exptime'][0], data['bg_flat_counts'][0], data['detxs'][0],
            data['t0'][0], data['t1'][0], data['mag_mcatbgsub_err_2'][0],
            data['flux'][0], data['mag_mcatbgsub_err_1'][0], data['flags'][0],
            data['mag_bgsub_err_2'][0], data['detrad'][0], data['cps'][0],
            data['flux_mcatbgsub_err'][0], data['flux_mcatbgsub'][0],
            mcat[band]['expt'][i],
            np.array(mcat[band]['ra'][i], dtype='float32'),
            np.array(mcat[band]['dec'][i], dtype='float32'),
            mcat[band][4]['mag'][i],
            mcat[band][4]['err'][i], mcat[band]['skybg'][i],
            mcat[band][7]['mag'][i], mcat[band][7]['err'][i])
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def datamaker(band, skypos, outfile, maglimit=20., margin=0.005,
              searchradius=0.1, radius=gt.aper2deg(4), annulus=[0.0083, 0.025],
              verbose=0):
    """
    Generate gAperture photometry for MCAT sources within a specified region.

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :param skypos: The right ascension and declination, in degrees.

    :type skypos: list

    :param outfile: Name of output file to make.

    :type outfile: str

    :param maglimit: Faint limit to use, in AB Mag.

    :type maglimit: float

    :param margin: The margin within which two sources are consider "the same,"
        in degrees.

    :type margin: float

    :param searchradius: The radius within which to search for sources, degrees.

    :type searchradius: float

    :param radius: The size of the aperture to measure fluxes with, in degrees.

    :type radius: float

    :param annulus: The inner and outer radii of the background annulus
        in degrees.

    :type annulus: float

    :param verbose: Verbosity level, a value of 0 is minimum verbosity.

    :type verbose: int
    """

    extant_objids = file_setup(outfile)

    if extant_objids == False:
        print('NOT RUNNING!!*!')
        return False

    uniques = dt.find_unique_sources(band, skypos[0], skypos[1], searchradius,
                                     maglimit=maglimit)

    if uniques is None:
        print('No sources at this position.')
        return

    for pos in uniques:
        mcat = dt.get_mcat_data(pos, margin)
        if not mcat:
            print('Nothing at {pos}.'.format(pos=pos))
            continue
        extant_objids = file_setup(outfile)
        for i, objid in enumerate(mcat['objid']):
            if mcat[band]['ra'][i] == -99. and mcat[band]['dec'][i] == -99.:
                print('No {b} source'.format(b=band))
                continue
            if objid in extant_objids:
                print('Already processed.')
                continue
            #exp = dt.exp_from_objid(objid)
            if mcat[band]['t0'][i] < 0:
                print('No MCAT exposure: skipping')
                continue
            print([mcat[band]['ra'][i], mcat[band]['dec'][i]])
            print([mcat[band]['t0'][i], mcat[band]['t1'][i]])
            data = gAperture(band, [mcat[band]['ra'][i], mcat[band]['dec'][i]],
                             radius, annulus=annulus, verbose=verbose,
                             coadd=True, trange=[mcat[band]['t0'][i],
                                                 mcat[band]['t1'][i]],
                             detsize=1.25)
            try:
                csv_construct = construct_row(i, band, objid, mcat, data)
                print(csv_construct)
                with open(outfile, 'ab') as csvfile:
                    spreadsheet = csv.writer(csvfile, delimiter=',',
                                             quotechar='|',
                                             quoting=csv.QUOTE_MINIMAL)
                    spreadsheet.writerow(csv_construct)
            except TypeError:
                continue

    return
# ------------------------------------------------------------------------------
