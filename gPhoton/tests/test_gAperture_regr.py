import unittest
import gPhoton.gMap as gm
import gPhoton.gFind as gf
import gPhoton.gAperture as ga

"""Regression tests for gAperture."""
class TestRegression(unittest.TestCase):
    def setUp(self):
        self.bands = ['NUV','FUV']
        self.skypos = [176.919525856024,0.255696872807351]
        self.tranges = [[ 765579107.995, 765579821.995 ],
                        [ 766525332.995, 766526576.995 ],
                        [ 919754986.995, 919755916.995 ]]
        self.trange = self.tranges[0]
        self.radius = 0.01
        self.annulus = [0.02,0.03]
        self.parser = ga.setup_parser()
        self.args = self.parser.parse_args()
        self.stepsz = 100.

    # FIXME: issue #95
    #def test_annulus_impact_on_mag(self):
    #    """Tests for a bug where the _not_ bg subtracted magnitude is
    #    different whether --annulus is specified or not."""
    #    for band in self.bands:
    #        out = ga.gAperture(band=band,skypos=self.skypos,
    #                                            radius=self.radius,coadd=True)
    #        out_ann = ga.gAperture(band=band,skypos=self.skypos,
    #                        radius=self.radius,annulus=self.annulus,coadd=True)
    #        self.assertAlmostEqual(out['mag'][0],out_ann['mag'][0])


    def test_basic_query_NUV(self):
        """Regtest the simplest valid query. (NUV)"""
        out = ga.gAperture(band='NUV',skypos=self.skypos,radius=self.radius)
        self.assertEqual(len(out['exptime']),6)
        # Regtest the exposure times
        self.assertAlmostEqual(out['exptime'][0],101.31698158)
        self.assertAlmostEqual(out['exptime'][1],635.41801554)
        self.assertAlmostEqual(out['exptime'][2],357.71295862)
        self.assertAlmostEqual(out['exptime'][3],1102.06018543)
        self.assertAlmostEqual(out['exptime'][4],96.51055404)
        self.assertAlmostEqual(out['exptime'][5],829.06433301)
        # Regtest the source magnitudes (no bg subtraction)
        self.assertAlmostEqual(out['mag'][0],17.79585028)
        self.assertAlmostEqual(out['mag'][1],17.78050662)
        self.assertAlmostEqual(out['mag'][2],17.83212313)
        self.assertAlmostEqual(out['mag'][3],13.16831294)
        self.assertAlmostEqual(out['mag'][4],17.83109209)
        self.assertAlmostEqual(out['mag'][5],17.79873666)

    def test_basic_query_FUV(self):
        """Regtest the simplest valid query. (FUV)"""
        out = ga.gAperture(band='FUV',skypos=self.skypos,radius=self.radius)
        self.assertEqual(len(out['exptime']),3)
        # Regtest the exposure times
        self.assertAlmostEqual(out['exptime'][0],1228.35811301)
        self.assertAlmostEqual(out['exptime'][1],107.05101862)
        self.assertAlmostEqual(out['exptime'][2],913.1313701)
        # Regtest the magnitudes (no bg subtraction)
        self.assertAlmostEqual(out['mag'][0],13.550352)
        self.assertAlmostEqual(out['mag'][1],19.16928608)
        self.assertAlmostEqual(out['mag'][2],19.06694976)

    def test_lcurve_query_NUV(self):
        """Regtest lightcurve over a specific time range. (NUV)"""
        out = ga.gAperture(band='NUV',skypos=self.skypos,radius=self.radius,
                           trange=self.tranges[1],stepsz=self.stepsz,
                           annulus=self.annulus)
        self.assertEqual(len(out['mag']),13)
        # Check that all of the bins have the right size (except the last one)
        for binsize in (out['t1']-out['t0'])[:-1]:
            self.assertAlmostEqual(binsize,self.stepsz)
        # Check the magnitudes (no bg subtraction)
        for i,mag in enumerate(
                               [ 17.19996981, 17.10094944, 16.99299101,
                                 16.86953449, 16.58609703, 15.64080472,
                                 12.8843604,  12.60970683, 12.63479485,
                                 11.95780728, 12.2635085,  12.553789,
                                 12.76446021 ]):
            self.assertAlmostEqual(out['mag'][i],mag)
        # Check the annulus background-subtracted fluxes (as magnitudes).
        for i,mag_bgsub in enumerate(
                                     [ 17.58702551, 17.44381185, 17.29197581,
                                       17.12539048, 16.78432337, 15.71999064,
                                       12.89102971, 12.61527998, 12.64022502,
                                       11.96186916, 12.268034,   12.55925879,
                                       12.77103454 ]):
            self.assertAlmostEqual(out['mag_bgsub'][i],mag_bgsub)
        # Check the annulus background values (counts, flat-weighted, and
        # within aperture).
        for i,bg_counts in enumerate([ 2110., 2089., 2051., 2007.,
                                       2073., 2089., 2293., 2461.,
                                       2341., 3248., 2735., 2530.,
                                       1057. ]):
            self.assertAlmostEqual(out['bg_counts'][i],bg_counts)
        for i,bg_flat_counts in enumerate(
                                          [ 1871.97722712, 1853.40782629,
                                            1821.46503959, 1780.35790353,
                                            1837.85420049, 1849.55875876,
                                            2032.06998278, 2184.33106164,
                                            2079.20667638, 2888.03488141,
                                            2429.30234821, 2245.24240439,
                                            941.65448052 ]):
            self.assertAlmostEqual(out['bg_flat_counts'][i],bg_flat_counts)
        for i,bg in enumerate(
                              [ 374.39544542, 370.68156526, 364.29300792,
                                356.07158071, 367.5708401,  369.91175175,
                                406.41399656, 436.86621233, 415.84133528,
                                577.60697628, 485.86046964, 449.04848088,
                                188.3308961 ]):
            self.assertAlmostEqual(out['bg'][i],bg)
        # Check the MCAT per-visit background subtracted fluxes.
        for i,mag_mcatbgsub in enumerate(
                                         [ 17.70991261, 17.55559111,
                                           17.39546401, 17.22086535,
                                           16.84606843, 15.74218895,
                                           12.8920306,  12.61565799,
                                           12.64088552, 11.96106793,
                                           12.26783162, 12.55944065,
                                           12.77132593 ]):
            self.assertAlmostEqual(out['mag_mcatbgsub'][i],mag_mcatbgsub)
        # Check the MCAT per-visit background values.
        for i,mcat_bg in enumerate(
                                   [ 467.93802676, 468.33843609, 468.76286039,
                                     468.86925093, 469.01876064, 468.8591898,
                                     467.19104717, 466.41687855, 466.27991533,
                                     463.84206572, 464.17653866, 463.93932338,
                                     196.65157325 ]):
            self.assertAlmostEqual(out['mcat_bg'][i],mcat_bg)

    def test_lcurve_query_FUV(self):
        """Regtest lightcurve over a specific time range. (FUV)"""
        out = ga.gAperture(band='FUV',skypos=self.skypos,radius=self.radius,
                           trange=self.tranges[1],stepsz=self.stepsz,
                           annulus=self.annulus)
        self.assertEqual(len(out['mag']),13)
        # Check that all of the bins have the right size (except the last one)
        for binsize in (out['t1']-out['t0'])[:-1]:
            self.assertAlmostEqual(binsize,self.stepsz)
        # Check the magnitudes (no bg subtraction)
        for i,mag in enumerate(
                               [ 18.87976244, 18.50287132, 18.57996131,
                                 18.49695605, 18.07555158, 16.86195021,
                                 13.17351786, 13.04236207, 13.0682284,
                                 12.09391959, 12.64711125, 13.23408717,
                                 13.54495825 ]):
            self.assertAlmostEqual(out['mag'][i],mag)
        # Check the annulus background-subtracted fluxes (as magnitudes).
        for i,mag_bgsub in enumerate(
                                     [ 19.22936765, 18.7231948,  18.84792398,
                                       18.73319054, 18.24969186, 16.92584474,
                                       13.17918486, 13.04696094, 13.07317331,
                                       12.10031926, 12.6511912, 13.23874903,
                                       13.55070637 ]):
            self.assertAlmostEqual(out['mag_bgsub'][i],mag_bgsub)
        # Check the annulus background values (counts, flat-weighted, and
        # within aperture).
        for i,bg_counts in enumerate([ 152., 143., 159., 153., 171., 201.,
                                       546., 501., 530., 1672., 644., 427.,
                                       170. ]):
            self.assertAlmostEqual(out['bg_counts'][i],bg_counts)
        for i,bg_flat_counts in enumerate(
                                          [ 127.47075602, 120.36695248,
                                            133.53371262, 128.87274145,
                                            144.00702825, 169.82327414,
                                            461.82457166, 423.16061006,
                                            444.22955368, 1406.82880576,
                                            540.05583777, 359.42448472,
                                            141.12829624 ]):
            self.assertAlmostEqual(out['bg_flat_counts'][i],bg_flat_counts)
        for i,bg in enumerate(
                              [ 25.4941512,   24.0733905,   26.70674252,
                                25.77454829,  28.80140565,  33.96465483,
                                92.36491433,  84.63212201,  88.84591074,
                                281.36576115, 108.01116755, 71.88489694,
                                28.22565925 ]):
            self.assertAlmostEqual(out['bg'][i],bg)
        # Check the MCAT per-visit background subtracted fluxes.
        for i,mag_mcatbgsub in enumerate(
                                         [ 19.65979013, 18.99112506,
                                           19.11455631, 18.98186517,
                                           18.37973659, 16.95236449,
                                           13.17642514, 13.04493815,
                                           13.07086667, 12.09499428,
                                           12.64890063, 13.23716149,
                                           13.5490537 ]):
            self.assertAlmostEqual(out['mag_mcatbgsub'][i],mag_mcatbgsub)
        # Check the MCAT per-visit background values.
        for i,mcat_bg in enumerate(
                                   [ 47.45846584, 47.4727517,  47.48119299,
                                     47.4810555,  47.48947085, 47.48532284,
                                     47.44526692, 47.45122304, 47.45240596,
                                     47.36531334, 47.42138539, 47.44004754,
                                     20.12567625 ]):
            self.assertAlmostEqual(out['mcat_bg'][i],mcat_bg)

    def test_coadd_query_NUV(self):
        """Regtest a coadded magnitude. (NUV)"""
        out = ga.gAperture(band='NUV',skypos=self.skypos,radius=self.radius,
                           tranges=self.tranges,coadd=True)
        # This query does *not* coadd, so we can compare the components that
        # went into the coadd values.
        out2 = ga.gAperture(band='NUV',skypos=self.skypos,radius=self.radius,
                           tranges=self.tranges,coadd=False)
        self.assertEqual(len(out['mag']), 1)
        self.assertEqual(len(out2['mag']), 3)
        # Exposure time of coadd must be sum of individual visits.
        self.assertEqual(sum(out2['exptime']), out['exptime'][0])
        # Check value of coadded magnitude.
        self.assertAlmostEqual(out['mag'][0],14.06593323)
        # Check that the uncertaintiy in the mag is less than the individual
        # visits.
        for x in out2['mag_err_1']:
            self.assertLess(out['mag_err_1'][0], x)
        for x in out2['mag_err_2']:
            self.assertLess(out['mag_err_2'][0], x)

    def test_coadd_query_FUV(self):
        """Regtest a coadded magnitude. (FUV)"""
        out = ga.gAperture(band='FUV',skypos=self.skypos,radius=self.radius,
                           tranges=self.tranges,coadd=True)
        # This query does *not* coadd, so we can compare the components that
        # went into the coadd values.
        out2 = ga.gAperture(band='FUV',skypos=self.skypos,radius=self.radius,
                           tranges=self.tranges,coadd=False)
        self.assertEqual(len(out['mag']), 1)
        self.assertEqual(len(out2['mag']), 2)
        # Exposure time of coadd must be sum of individual visits.
        self.assertEqual(sum(out2['exptime']), out['exptime'][0])
        # Check value of coadded magnitude.
        self.assertAlmostEqual(out['mag'][0],14.14882565)
        # Check that the uncertaintiy in the mag is less than the individual
        # visits.
        for x in out2['mag_err_1']:
            self.assertLess(out['mag_err_1'][0], x)
        for x in out2['mag_err_2']:
            self.assertLess(out['mag_err_2'][0], x)

suite = unittest.TestLoader().loadTestsFromTestCase(TestRegression)
unittest.TextTestRunner(verbosity=2).run(suite)
