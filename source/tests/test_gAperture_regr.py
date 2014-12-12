import unittest
import gMap as gm
import gFind as gf
import gAperture as ga

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

    def test_basic_query_NUV(self):
        """Regtest the simplest valid query. (NUV)"""
        out = ga.gAperture(band='NUV',skypos=self.skypos,radius=self.radius)
        self.assertEqual(len(out['exptime']),4)
        # Regtest the exposure times
        self.assertAlmostEqual(out['exptime'][0],634.2852623)
        self.assertAlmostEqual(out['exptime'][1],1099.98059419)
        self.assertAlmostEqual(out['exptime'][2],457.07447477)
        self.assertAlmostEqual(out['exptime'][3],370.61241714)
        # Regtest the source magnitudes (no bg subtraction)
        self.assertAlmostEqual(out['mag'][0],17.77856936)
        self.assertAlmostEqual(out['mag'][1],13.16626222)
        self.assertAlmostEqual(out['mag'][2],17.821252)
        self.assertAlmostEqual(out['mag'][3],17.76766803)

    def test_basic_query_FUV(self):
        """Regtest the simplest valid query. (FUV)"""
        out = ga.gAperture(band='FUV',skypos=self.skypos,radius=self.radius)
        self.assertEqual(len(out['exptime']),2)
        # Regtest the exposure times
        self.assertAlmostEqual(out['exptime'][0],1227.36996722)
        self.assertAlmostEqual(out['exptime'][1],406.73730805)
        # Regtest the magnitudes (no bg subtraction)
        self.assertAlmostEqual(out['mag'][0],13.54947823)
        self.assertAlmostEqual(out['mag'][1],19.03458803)

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
        for i,mag in enumerate([ 17.19836281,  17.09897678,  17.00228114,
                                 16.86768186,  16.58436345,  15.63889791,
                                 12.88250218,  12.61903232,  12.64400935,
                                 11.96701482,  12.27271088,  12.56295656,
                                 12.7914485 ]):
            self.assertAlmostEqual(out['mag'][i],mag)
        # w/ background subtraction: What this might really be testing is the
        #  efficiency of the Monte Carlo annulus area estimate...
        for i,bgsubmag in enumerate([ 17.58541851,  17.44183919,  17.30126593,
                                      17.12353785,  16.78258979,  15.71808383,
                                      12.88917149,  12.62460547,  12.64943953,
                                      11.9710767 ,  12.27723637,  12.56842635,
                                      12.79805813]):
            self.assertAlmostEqual(out['mag_bgsub_cheese'][i],bgsubmag)

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
        for i,mag in enumerate([ 18.87932679,  18.48008581,  18.55708448,
                                 18.47418146,  18.07477455,  16.8611657 ,
                                 13.17292039,  13.05277282,  13.06770171,
                                 12.10430201,  12.65757588,  13.23356214,
                                 13.55337665]):
            self.assertAlmostEqual(out['mag'][i],mag)
        # w/ background subtraction: What this might really be testing is the
        #  efficiency of the Monte Carlo annulus area estimate...
        for i,bgsubmag in enumerate([ 19.22893201,  18.70040929,  18.82504715,
                                      18.71041594,  18.24891483,  16.92506023,
                                      13.1785874 ,  13.05737169,  13.07264662,
                                      12.11070168,  12.66165582,  13.238224  ,
                                      13.55910122]):
            self.assertAlmostEqual(out['mag_bgsub_cheese'][i],bgsubmag)

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

    def test_coadd_query_NUV(self):
        """Regtest a coadded magnitude. (NUV)"""
        out = ga.gAperture(band='NUV',skypos=self.skypos,radius=self.radius,
                           tranges=self.tranges,coadd=True)
        self.assertAlmostEqual(out['mag'][0],14.064324896330707)

    def test_coadd_query_FUV(self):
        """Regtest a coadded magnitude. (FUV)"""
        out = ga.gAperture(band='FUV',skypos=self.skypos,radius=self.radius,
                           tranges=self.tranges,coadd=True)
        self.assertAlmostEqual(out['mag'][0],13.855242986154447)

suite = unittest.TestLoader().loadTestsFromTestCase(TestRegression)
unittest.TextTestRunner(verbosity=2).run(suite)
