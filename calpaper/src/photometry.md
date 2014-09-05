##Photometry
The canonical GALEX pipeline used the Source Extractor (SExtractor) toolkit to extract photometric information from images. gAperture uses naive aperture photometry. These methods needs to be compared. The most straightforward way to do this is to compare gAperture flux measurements to their corresponding MCAT or GCAT values for a large number of sources across the sky. Do this efficiently will require enhancement of gAperture with the capability to determine and automatically use "optimum" aperture parameters based upon the MCAT/GCAT entry for the target in question.

A list of potential White Dwarf Standards appears in Table 5 of Morrissey, Patrick, et al. "The calibration and data products of GALEX." The Astrophysical Journal Supplement Series 173.2 (2007): 682.

Morrissey, et al. suggests a 12" (0.003 degree) aperture as a good generic value.

"For both detectors, the onset of saturation is at approximately mAB ~ 15 (corresponding to 34 counts s^-1 FUV and 108 counts s^-1 NUV) _in the lowest-gain regions of the detector_."

###LDS749B
This is the primary white dwarf calibration standard for GALEX. It has nominal AB magnitudes of *14.71* in NUV and *15.57* in FUV. Note that this is past the saturation point in NUV, which needs to be accounted for in the analysis.

    ./gAperture.py --skypos [323.06766667,0.25400000] -a 0.005 --annulus [0.01,0.02] -b 'NUV' -v 2 -f 'LDS749B_NUV_rr.csv' --maxgap 1600 --minexp 1000 --hrbg --response --overwrite

In [222]: counts2mag(mag2counts(data['mag']-data['apcorr2'],'NUV').mean(),'NUV')Out[222]: 14.713075354913997

    ./gAperture.py --skypos [323.06766667,0.25400000] -a 0.005 --annulus [0.01,0.02] -b 'FUV' -v 2 -f 'LDS749B_FUV_rr.csv' --maxgap 1600 --minexp 1000 --hrbg --response --overwrite

In [240]: counts2mag(mag2counts(data['mag']-data['apcorr2'],'FUV').mean(),'FUV')
Out[240]: 15.414742366262452

    import gAperture
    skypos = [323.06766667,0.25400000]
    ra0,dec0 = skypos
    radius = 0.005
    annulus = [0.01,0.02]
    band = 'FUV'

    data = gAperture.gAperture('FUV',[323.06766667,0.25400000],0.005,
                               annulus=[0.01,0.02],stepsz=1500,verbose=2)
    data = gAperture.gAperture('NUV',[323.06766667,0.25400000],0.005,
                               annulus=[0.01,0.02],stepsz=1500,verbose=2,
                               trange=[0,870103402.995])

    %pylab
    import galextools as gt
    catmag,band = 15.57,'FUV'
    ap2 = gt.apcorrect2(band,radius)
    import gphoton_utils as gu
    plt.figure()
    plt.gca().invert_yaxis()
    plt.plot(data['exptime'],data['mag']-ap2,'.')
    err = gu.model_errors(15.57,'FUV',sigma=5)
    plt.plot(np.arange(1,1600),err[0])
    plt.plot(np.arange(1,1600),err[1])
    plt.plot(np.arange(0,1600),np.ones(1600)*catmag)

    data['dist']=sqrt((data['detxs']-400.)**2. + (data['detys']-400.)**2.)

    import gQuery as gq
    query = 'http://masttest.stsci.edu/portal/Mashup/MashupQuery.asmx/GalexPhotonListQueryTest?query=select ra, dec, nuv_mag, fuv_mag, fov_radius, nuv_skybg, fuv_skybg, nuv_fwhm_world, fuv_fwhm_world, vpe.fexptime, vpe.nexptime, n_zpmag, f_zpmag, nuv_mag_aper_4, fuv_mag_aper_4 from Gr6plus7.Dbo.visitphotoobjall as vpo inner join Gr6plus7.Dbo.visitphotoextract as vpe on vpo.photoextractid=vpe.photoextractid inner join gr6plus7.dbo.fGetNearbyVisitObjEq(323.06766667, 0.254, 0.1) as nb on vpo.objid=nb.objid&format=json&timeout={}'
    out = np.array(gq.getArray(query))
    LDS749B = {'ra':out[:,0],'dec':out[:,1],'nmag':out[:,2],'fmag':out[:,3],'fexpt':out[:,9],'nexpt':out[:,10],'nzp':out[:,11],'fzp':out[:,12],'n4':out[:,13],'f4':out[:,14]}
    ix=np.where(LDS749B['fmag']>0)
    # The FUV catalog magnitude is converted to the true magnitude using the
    # zero-point conversion (see header card
    # "F_ZPMAG": fuv_mag = FUV_MAG_AUTO + F_ZPMAG)
    plt.plot(LDS749B['fexpt'][ix],LDS749B['fmag'][ix],'x')
    plt.xlabel('Exposure time (s)')
    plt.ylabel('AB Mag')
    plt.title('LDS749B - FUV')


-----------

    plt.figure()
    catmag,band = 14.71,'NUV'
    err = gu.model_errors(catmag,'FUV',sigma=5)
    plt.plot(np.arange(1,1600),err[0])
    plt.plot(np.arange(1,1600),err[1])
    plt.plot(np.arange(0,1600),np.ones(1600)*catmag)
    ix=np.where(LDS749B['nmag']>0)
    plt.plot(LDS749B['nexpt'][ix],LDS749B['nmag'][ix],'.')
    plt.xlabel('Exposure time (s)')
    plt.ylabel('AB Mag')
    plt.title('LDS749B - NUV')


