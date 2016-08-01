import unittest
import gPhoton.gMap as gm
import gPhoton.gFind as gf
import gPhoton.gAperture as ga

from gPhoton.gAperture import setup_parser as ga_setup_parser

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
        self.parser = ga_setup_parser()
        self.args = self.parser.parse_args()
        self.stepsz = 100.

    # FIXME: issue #95
    #def test_annulus_impact_on_mag(self):
    #    """Tests for a bug where the _not_ bg subtracted magnitude is
    #    different whether --annulus is specified or not."""
    #    for band in self.bands:
    #        out = ga(band=band,skypos=self.skypos,
    #                                            radius=self.radius,coadd=True)
    #        out_ann = ga(band=band,skypos=self.skypos,
    #                        radius=self.radius,annulus=self.annulus,coadd=True)
    #        self.assertAlmostEqual(out['mag'][0],out_ann['mag'][0])


    def test_basic_query_NUV(self):
        """Regtest the simplest valid query. (NUV)"""
        out = ga(band='NUV',skypos=self.skypos,radius=self.radius)
        self.assertEqual(len(out['exptime']),6)
        # Regtest the exposure times
        for i, expt in enumerate([ 102.20581307, 636.30367737, 358.59845156,
                                   1102.93999215, 97.40319135, 829.95461701]):
            self.assertAlmostEqual(out['exptime'][i],expt)
        # Regtest the source magnitudes (no bg subtraction)
        for i, mag in enumerate([ 17.80533367, 17.78201889, 17.83480747,
                                  13.16917937, 17.84108803, 17.79990195]):
            self.assertAlmostEqual(out['mag'][i],mag)

    def test_basic_query_FUV(self):
        """Regtest the simplest valid query. (FUV)"""
        out = ga(band='FUV',skypos=self.skypos,radius=self.radius)
        self.assertEqual(len(out['exptime']),3)
        # Regtest the exposure times
        for i, expt in enumerate([ 1229.34564077, 108.04223136, 914.1187588 ]):
            self.assertAlmostEqual(out['exptime'][i],expt)
        # Regtest the magnitudes (no bg subtraction)
        for i, mag in enumerate([ 13.55122452,  19.17929293,  19.07231177]):
            self.assertAlmostEqual(out['mag'][i],mag)

    def test_lcurve_query_NUV(self):
        """Regtest lightcurve over a specific time range. (NUV)"""
        out = ga(band='NUV',skypos=self.skypos,radius=self.radius,
                           trange=self.tranges[1],stepsz=self.stepsz,
                           annulus=self.annulus)
        self.assertEqual(len(out['mag']),13)
        # Check that all of the bins have the right size (except the last one)
        for binsize in (out['t1']-out['t0'])[:-1]:
            self.assertAlmostEqual(binsize,self.stepsz)
        # Check the magnitudes (no bg subtraction)
        for i,mag in enumerate([ 17.21087962,  17.1118642 ,  17.00391571,
                                 16.88044162,  16.59700633,  15.65171105,
                                 12.89521826,  12.62064997,  12.64563419,
                                 11.96873523,  12.27441822,  12.56467381,
                                 12.79256972]):
            self.assertAlmostEqual(out['mag'][i],mag)
        # Check the annulus background-subtracted fluxes (as magnitudes).
        for i,mag_bgsub in enumerate([ 17.59793532,  17.4547266 ,  17.30290051,
                                       17.13629761,  16.79523266,  15.73089697,
                                       12.90188757,  12.62622312,  12.65106437,
                                       11.9727971 ,  12.27894371,  12.5701436 ,
                                       12.79917936]):
            self.assertAlmostEqual(out['mag_bgsub'][i],mag_bgsub)
        # Check the annulus background values (counts, flat-weighted, and
        # within aperture).
        for i,bg_counts in enumerate([ 2110., 2089., 2051., 2007.,
                                       2073., 2089., 2293., 2461.,
                                       2341., 3248., 2735., 2530.,
                                       1085. ]):
            self.assertAlmostEqual(out['bg_counts'][i],bg_counts)
        for i,bg_flat_counts in enumerate([ 1871.97722712,  1853.40782629,
                                            1821.46503959,  1780.35790353,
                                            1837.85420049,  1849.55875876,
                                            2032.06998278,  2184.33106164,
                                            2079.20667638,  2888.03488141,
                                            2429.30234821,  2245.24240439,
                                            966.42340298]):
            self.assertAlmostEqual(out['bg_flat_counts'][i],bg_flat_counts)
        for i,bg in enumerate([ 374.39544542,  370.68156526,  364.29300792,
                                356.07158071,  367.5708401 ,  369.91175175,
                                406.41399656,  436.86621233,  415.84133528,
                                577.60697628,  485.86046964,  449.04848088,
                                193.2846806 ]):
            self.assertAlmostEqual(out['bg'][i],bg)
        # Check the MCAT per-visit background subtracted fluxes.
        for i,mag_mcatbgsub in enumerate([ 17.72741551, 17.57222571,
            17.41132697, 17.23596895, 16.85994807, 15.75416867, 12.90296583,
            12.62666159, 12.65178614,  11.97202891,  12.27878508,  12.57038255,
            12.7996161 ]):
            self.assertAlmostEqual(out['mag_mcatbgsub'][i],mag_mcatbgsub)
        # Check the MCAT per-visit background values.
        for i,mcat_bg in enumerate([ 472.66371408,  473.07032179,  473.50336992,
            473.60317282, 473.75513792,  473.59266252,  471.88661096,
            471.14166756, 470.95829134,  468.53421591,  468.86419727,
            468.61383664, 206.01485305]):
            self.assertAlmostEqual(out['mcat_bg'][i],mcat_bg)

    def test_lcurve_query_FUV(self):
        """Regtest lightcurve over a specific time range. (FUV)"""
        out = ga(band='FUV',skypos=self.skypos,radius=self.radius,
                           trange=self.tranges[1],stepsz=self.stepsz,
                           annulus=self.annulus)
        self.assertEqual(len(out['mag']),13)
        # Check that all of the bins have the right size (except the last one)
        for binsize in (out['t1']-out['t0'])[:-1]:
            self.assertAlmostEqual(binsize,self.stepsz)
        # Check the magnitudes (no bg subtraction)
        for i,mag in enumerate([ 18.89068118,  18.51378867,  18.59087309,
            18.50787086, 18.08646446, 16.87285734, 13.1844053, 13.05328141,
            13.07911962,  12.10483921,  12.65802432,  13.24499784,
            13.57836495]):
            self.assertAlmostEqual(out['mag'][i],mag)
        # Check the annulus background-subtracted fluxes (as magnitudes).
        for i,mag_bgsub in enumerate([ 19.24028639,  18.73411215,  18.85883575,
            18.74410534, 18.26060474,  16.93675188,  13.1900723 ,  13.05788028,
            13.08406453,  12.11123888,  12.66210426,  13.2496597 ,
            13.58408952]):
            self.assertAlmostEqual(out['mag_bgsub'][i],mag_bgsub)
        # Check the annulus background values (counts, flat-weighted, and
        # within aperture).
        for i,bg_counts in enumerate([  152.,   143.,   159.,   153.,   171.,
            201.,   546.,   501., 530.,  1672.,   644.,   427.,   172.]):
            self.assertAlmostEqual(out['bg_counts'][i],bg_counts)
        for i,bg_flat_counts in enumerate([  127.47075602,   120.36695248,
            133.53371262,   128.87274145, 144.00702825,   169.82327414,
            461.82457166,   423.16061006, 444.22955368,  1406.82880576,
            540.05583777,   359.42448472, 142.78094526]):
            self.assertAlmostEqual(out['bg_flat_counts'][i],bg_flat_counts)
        for i,bg in enumerate([  25.4941512 ,   24.0733905 ,   26.70674252,
            25.77454829, 28.80140565,   33.96465483,   92.36491433,
            84.63212201, 88.84591074,  281.36576115,  108.01116755,
            71.88489694, 28.55618905]):
            self.assertAlmostEqual(out['bg'][i],bg)
        # Check the MCAT per-visit background subtracted fluxes.
        for i,mag_mcatbgsub in enumerate([ 19.68230643,  19.00829092,
            19.13246775,  18.99897381, 18.39420172,  16.964224  ,  13.18734192,
            13.05588356, 13.08178451,  12.10592476,  12.65983179,  13.24810326,
            13.58258862]):
            self.assertAlmostEqual(out['mag_mcatbgsub'][i],mag_mcatbgsub)
        # Check the MCAT per-visit background values.
        for i,mcat_bg in enumerate([  47.93814125, 47.95251012, 47.96079057,
            47.96078538, 47.96920073, 47.96475699, 47.92342736, 47.93085175,
            47.93080506, 47.84408598, 47.90043561, 47.91918077, 21.08372221]):
            self.assertAlmostEqual(out['mcat_bg'][i],mcat_bg)

    def test_coadd_query_NUV(self):
        """Regtest a coadded magnitude. (NUV)"""
        out = ga(band='NUV',skypos=self.skypos,radius=self.radius,
                           tranges=self.tranges,coadd=True)
        # This query does *not* coadd, so we can compare the components that
        # went into the coadd values.
        out2 = ga(band='NUV',skypos=self.skypos,radius=self.radius,
                           tranges=self.tranges,coadd=False)
        self.assertEqual(len(out['mag']), 1)
        self.assertEqual(len(out2['mag']), 3)
        # Exposure time of coadd must be sum of individual visits.
        self.assertEqual(sum(out2['exptime']), out['exptime'][0])
        # Check value of coadded magnitude.
        self.assertAlmostEqual(out['mag'][0],14.067056127080908)
        # Check that the uncertaintiy in the mag is less than the individual
        # visits.
        for x in out2['mag_err_1']:
            self.assertLess(out['mag_err_1'][0], x)
        for x in out2['mag_err_2']:
            self.assertLess(out['mag_err_2'][0], x)

    def test_coadd_query_FUV(self):
        """Regtest a coadded magnitude. (FUV)"""
        out = ga(band='FUV',skypos=self.skypos,radius=self.radius,
                           tranges=self.tranges,coadd=True)
        # This query does *not* coadd, so we can compare the components that
        # went into the coadd values.
        out2 = ga(band='FUV',skypos=self.skypos,radius=self.radius,
                           tranges=self.tranges,coadd=False)
        self.assertEqual(len(out['mag']), 1)
        self.assertEqual(len(out2['mag']), 2)
        # Exposure time of coadd must be sum of individual visits.
        self.assertEqual(sum(out2['exptime']), out['exptime'][0])
        # Check value of coadded magnitude.
        self.assertAlmostEqual(out['mag'][0],14.1498456983)
        # Check that the uncertaintiy in the mag is less than the individual
        # visits.
        for x in out2['mag_err_1']:
            self.assertLess(out['mag_err_1'][0], x)
        for x in out2['mag_err_2']:
            self.assertLess(out['mag_err_2'][0], x)

suite = unittest.TestLoader().loadTestsFromTestCase(TestRegression)
unittest.TextTestRunner(verbosity=2).run(suite)
