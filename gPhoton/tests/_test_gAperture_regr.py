from __future__ import absolute_import, division, print_function
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

    def test_annulus_impact_on_mag(self):
        """Tests for a bug where the _not_ bg subtracted magnitude is
        different whether --annulus is specified or not."""
        for band in self.bands:
            out = ga(band=band,skypos=self.skypos,
                                                radius=self.radius,coadd=True)
            out_ann = ga(band=band,skypos=self.skypos,
                            radius=self.radius,annulus=self.annulus,coadd=True)
            self.assertAlmostEqual(out['mag'][0],out_ann['mag'][0])

    def test_basic_query_NUV(self):
        """Regtest the simplest valid query. (NUV)"""
        out = ga(band='NUV',skypos=self.skypos,radius=self.radius)
        self.assertEqual(len(out['exptime']),6)
        # Regtest the exposure times
        for i, expt in enumerate([ 102.20581307, 636.30367737, 358.59845156,
                                   1102.93999215, 97.40319135, 829.95018001]):
            self.assertAlmostEqual(out['exptime'][i],expt)
        # Regtest the source magnitudes (no bg subtraction)
        for i, mag in enumerate([ 17.80533367, 17.78201889, 17.83480747,
                                  13.16917937, 17.84108803, 17.79989614]):
            self.assertAlmostEqual(out['mag'][i],mag)

    def test_basic_query_FUV(self):
        """Regtest the simplest valid query. (FUV)"""
        out = ga(band='FUV',skypos=self.skypos,radius=self.radius)
        self.assertEqual(len(out['exptime']),3)
        # Regtest the exposure times
        for i, expt in enumerate([ 1229.34564077, 108.04223136, 914.11382299 ]):
            self.assertAlmostEqual(out['exptime'][i],expt)
        # Regtest the magnitudes (no bg subtraction)
        for i, mag in enumerate([ 13.55122452,  19.17929293,  19.07230591]):
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
        for i,mag in enumerate([17.210825  , 17.11186399, 17.00391599,
            16.88044217, 16.59700688, 15.65171036, 12.89538741, 12.6205782 ,
            12.64566182, 11.96872796, 12.27440602, 12.56467384, 12.79249614]):
            self.assertAlmostEqual(out['mag'][i],mag)
        # Check the annulus background-subtracted fluxes (as magnitudes).
        for i,mag_bgsub in enumerate([17.59766376, 17.45491639, 17.30290078,
            17.13615794, 16.79533561, 15.73089628, 12.90205776, 12.62615098,
            12.65109214, 11.9727898 , 12.27893147, 12.57014363, 12.79910532]):
            self.assertAlmostEqual(out['mag_bgsub'][i],mag_bgsub)
        # Check the annulus background values (counts, flat-weighted, and
        # within aperture).
        for i,bg_counts in enumerate([2109., 2090., 2051., 2006., 2074.,
            2089., 2293., 2461., 2341., 3248., 2735., 2530., 1085.]):
            self.assertAlmostEqual(out['bg_counts'][i],bg_counts)
        for i,bg_flat_counts in enumerate([1871.10387552, 1854.28117788,
            1821.46503959, 1779.49259736, 1838.71950666, 1849.55875876,
            2032.06998278, 2184.33106164, 2079.20667638, 2888.03488141,
            2429.30234821, 2245.24240439, 966.42340298]):
            self.assertAlmostEqual(out['bg_flat_counts'][i],bg_flat_counts)
        for i,bg in enumerate([374.2207751 , 370.85623558, 364.29300792,
            355.89851947, 367.74390133, 369.91175175, 406.41399656,
            436.86621233, 415.84133528, 577.60697628, 485.86046964,
            449.04848088, 193.2846806 ]):
            self.assertAlmostEqual(out['bg'][i],bg)
        # Check the MCAT per-visit background subtracted fluxes.
        for i,mag_mcatbgsub in enumerate([17.72732761, 17.57222539, 17.41132737,
            17.23596971, 16.85994877, 15.75416792, 12.90313619, 12.62658942,
            12.65181393, 11.97202162, 12.27877284, 12.57038258, 12.79954204]):
            self.assertAlmostEqual(out['mag_mcatbgsub'][i],mag_mcatbgsub)
        # Check the MCAT per-visit background values.
        for i,mcat_bg in enumerate([472.63993535, 473.07023196, 473.5034897 ,
            473.60341237, 473.75537747, 473.59236308, 471.88655108, 471.14208677,
            470.95781224, 468.53454529, 468.86386789, 468.61389653, 206.01527226]):
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
        for i,mag in enumerate([18.8906269 , 18.51378855, 18.59087303, 18.50787103, 18.08646446,
            16.87285722, 13.18448178, 13.05321392, 13.07925738, 12.10481064,
            12.65807191, 13.24499874, 13.57811565]):
            self.assertAlmostEqual(out['mag'][i],mag)
        # Check the annulus background-subtracted fluxes (as magnitudes).
        for i,mag_bgsub in enumerate([19.24023211, 18.73411203, 18.8569994 , 18.74576033, 18.26060474,
            16.93675176, 13.19014918, 13.05781251, 13.08420293, 12.11121014,
            12.66215203, 13.2496606 , 13.58380449]):
            self.assertAlmostEqual(out['mag_bgsub'][i],mag_bgsub)
        # Check the annulus background values (counts, flat-weighted, and
        # within aperture).
        for i,bg_counts in enumerate([ 152.,  143.,  158.,  154.,  171.,
            201.,  546.,  501.,  530., 1672.,  644.,  427.,  171.]):
            self.assertAlmostEqual(out['bg_counts'][i],bg_counts)
        for i,bg_flat_counts in enumerate([ 127.47075602,  120.36695248,  132.72623196,  129.68022211,
            144.00702825,  169.82327414,  461.82457166,  423.16061006,
            444.22955368, 1406.82880576,  540.05583777,  359.42448472,
            141.92463986]):
            self.assertAlmostEqual(out['bg_flat_counts'][i],bg_flat_counts)
        for i,bg in enumerate([ 25.4941512 ,  24.0733905 ,  26.54524639,  25.93604442,
            28.80140565,  33.96465483,  92.36491433,  84.63212201,
            88.84591074, 281.36576115, 108.01116755,  71.88489694,
            28.38492797]):
            self.assertAlmostEqual(out['bg'][i],bg)
        # Check the MCAT per-visit background subtracted fluxes.
        for i,mag_mcatbgsub in enumerate([19.68219389, 19.00829074, 19.13246765, 18.99897409, 18.39420172,
            16.96422387, 13.18741861, 13.05581591, 13.08192262, 12.10589617,
            12.65987946, 13.24810416, 13.58233835]):
            self.assertAlmostEqual(out['mag_mcatbgsub'][i],mag_mcatbgsub)
        # Check the MCAT per-visit background values.
        for i,mcat_bg in enumerate([47.9357446 , 47.95250493, 47.96078798, 47.96079316, 47.96920073,
            47.9647518 , 47.92345589, 47.93083878, 47.9307843 , 47.84410414,
            47.90043042, 47.91917558, 21.08374297]):
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
