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
        self.assertAlmostEqual(out['exptime'][0],102.1226448)
        self.assertAlmostEqual(out['exptime'][1],635.79861454)
        self.assertAlmostEqual(out['exptime'][2],358.3098205)
        self.assertAlmostEqual(out['exptime'][3],1102.0230354)
        self.assertAlmostEqual(out['exptime'][4],97.32780661)
        self.assertAlmostEqual(out['exptime'][5],829.32579227)
        # Regtest the source magnitudes (no bg subtraction)
        self.assertAlmostEqual(out['mag'][0],17.80444981)
        self.assertAlmostEqual(out['mag'][1],17.78115675)
        self.assertAlmostEqual(out['mag'][2],17.83393323)
        self.assertAlmostEqual(out['mag'][3],13.16827634)
        self.assertAlmostEqual(out['mag'][4],17.84024741)
        self.assertAlmostEqual(out['mag'][5],17.79907901)

    def test_basic_query_FUV(self):
        """Regtest the simplest valid query. (FUV)"""
        out = ga.gAperture(band='FUV',skypos=self.skypos,radius=self.radius)
        self.assertEqual(len(out['exptime']),3)
        # Regtest the exposure times
        self.assertAlmostEqual(out['exptime'][0],1229.27537977)
        self.assertAlmostEqual(out['exptime'][1],108.03763929)
        self.assertAlmostEqual(out['exptime'][2],917.04785491)
        # Regtest the magnitudes (no bg subtraction)
        self.assertAlmostEqual(out['mag'][0],13.55116246)
        self.assertAlmostEqual(out['mag'][1],19.17924678)
        self.assertAlmostEqual(out['mag'][2],19.0715966)

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
        for i,mag in enumerate([ 17.21694316, 17.11307157, 17.00687142,
                                 16.88015838, 16.59671077, 15.66001062,
                                 12.90278139, 12.61704438, 12.64491365,
                                 11.96610893, 12.27605849, 12.57045386,
                                 12.76668811 ]):
            self.assertAlmostEqual(out['mag'][i],mag)
        # Check the annulus background-subtracted fluxes (as magnitudes).
        for i,mag_bgsub in enumerate([ 17.60929357, 17.45840616, 17.30776929,
                                       17.13685966, 16.79604557, 15.7404064,
                                       12.90954336, 12.62263757, 12.6503745,
                                       11.97017948, 12.28061957, 12.57599662,
                                       12.7731486 ]):
            self.assertAlmostEqual(out['mag_bgsub'][i],mag_bgsub)
        # Check the annulus background values (counts, flat-weighted, and
        # within aperture).
        for i,bg_counts in enumerate([ 2110., 2089., 2051., 2007.,
                                       2073., 2089., 2293., 2461.,
                                       2341., 3248., 2735., 2530.,
                                       1057. ]):
            self.assertAlmostEqual(out['bg_counts'][i],bg_counts)
        for i,bg_flat_counts in enumerate([ 1881.16239462, 1861.18207487,
                                            1825.13672239, 1784.60944618,
                                            1846.24637584, 1861.00296681,
                                            2044.22407791, 2197.61381261,
                                            2090.54813136, 2898.65336221,
                                            2442.5370132,  2261.0528663,
                                            944.66126792 ]):
            self.assertAlmostEqual(out['bg_flat_counts'][i],bg_flat_counts)
        for i,bg in enumerate([ 376.23247892, 372.23641497, 365.02734448,
                                356.92188924, 369.24927517, 372.20059336,
                                408.84481558, 439.52276252, 418.10962627,
                                579.73067244, 488.50740264, 452.21057326,
                                188.93225358 ]):
            self.assertAlmostEqual(out['bg'][i],bg)
        # Check the MCAT per-visit background subtracted fluxes.
        for i,mag_mcatbgsub in enumerate([ 17.73718973, 17.57407121,
                                           17.4156312, 17.23557599,
                                           16.85957154, 15.76329305,
                                           12.91058331, 12.62303601,
                                           12.65106151, 11.96939464,
                                           12.28043197, 12.57619315,
                                           12.77356797]):
            self.assertAlmostEqual(out['mag_mcatbgsub'][i],mag_mcatbgsub)
        # Check the MCAT per-visit background values.
        for i,mcat_bg in enumerate([ 472.27884811, 472.68809896,
                                     473.12396211, 473.22441377,
                                     473.37736672, 473.21383515,
                                     471.49669346, 470.74690756,
                                     470.56233932, 468.12250625,
                                     468.45463265, 468.20264455,
                                     201.15752093]):
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
        for i,mag in enumerate([ 18.8824072, 18.51488516, 18.58549817,
                                 18.49929593, 18.08761381, 16.85940851,
                                 13.17413591, 13.066734, 13.08552948,
                                 12.0927487, 12.67582645, 13.21586764,
                                 13.61488068 ]):
            self.assertAlmostEqual(out['mag'][i],mag)
        # Check the annulus background-subtracted fluxes (as magnitudes).
        for i,mag_bgsub in enumerate([ 19.2232331, 18.73177926, 18.84716332,
                                       18.72850653, 18.25984486, 16.92103755,
                                       13.17962333, 13.0712827, 13.09040537,
                                       12.09893898, 12.67990986, 13.22032368,
                                       13.62084791 ]):
            self.assertAlmostEqual(out['mag_bgsub'][i],mag_bgsub)
        # Check the annulus background values (counts, flat-weighted, and
        # within aperture).
        for i,bg_counts in enumerate([ 152., 143., 159., 153., 171., 201.,
                                       546., 501., 530., 1672., 644., 427.,
                                       170. ]):
            self.assertAlmostEqual(out['bg_counts'][i],bg_counts)
        for i,bg_flat_counts in enumerate([ 125.69362214, 118.54809357,
                                            131.40057948, 126.41933067,
                                            142.39169606, 166.00599403,
                                            451.45002048, 413.37614878,
                                            435.43888582, 1376.0788091,
                                            531.69121736, 352.91165538,
                                            140.61735056 ]):
            self.assertAlmostEqual(out['bg_flat_counts'][i],bg_flat_counts)
        for i,bg in enumerate([ 25.13872443, 23.70961871, 26.2801159,
                                25.28386613, 28.47833921, 33.20119881,
                                90.2900041, 82.67522976, 87.08777716,
                                275.21576182, 106.33824347, 70.58233108,
                                28.12347011 ]):
            self.assertAlmostEqual(out['bg'][i],bg)
        # Check the MCAT per-visit background subtracted fluxes.
        for i,mag_mcatbgsub in enumerate([ 19.66522207, 19.01002047,
                                           19.1236305, 18.98552466,
                                           18.39572795, 16.94960235,
                                           13.17704485, 13.06936863,
                                           13.08821018, 12.09382223,
                                           12.67766382, 13.21889073,
                                           13.61924911 ]):
            self.assertAlmostEqual(out['mag_mcatbgsub'][i],mag_mcatbgsub)
        # Check the MCAT per-visit background values.
        for i,mcat_bg in enumerate([ 47.93542943, 47.94986719, 47.95818734,
                                     47.95818213, 47.96663783, 47.96217278,
                                     47.92064499, 47.92810498, 47.92805806,
                                     47.84092321, 47.89754301, 47.91637805,
                                     20.60352662]):
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
        self.assertAlmostEqual(out['mag'][0],14.24326228)
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
        self.assertAlmostEqual(out['mag'][0],14.20268608)
        # Check that the uncertaintiy in the mag is less than the individual
        # visits.
        for x in out2['mag_err_1']:
            self.assertLess(out['mag_err_1'][0], x)
        for x in out2['mag_err_2']:
            self.assertLess(out['mag_err_2'][0], x)

suite = unittest.TestLoader().loadTestsFromTestCase(TestRegression)
unittest.TextTestRunner(verbosity=2).run(suite)
