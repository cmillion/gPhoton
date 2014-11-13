import unittest
import gMap as gm
import gFind as gf
import gAperture as ga

"""Unit tests for the command line parameters and checking for the gPhoton
command line utilities [gPhoton.py, gAperture.py,gMap.py, gFind.py].
"""
class TestgApertureArguments(unittest.TestCase):
    def setUp(self):
        self.bands = ['NUV','FUV']
        self.skypos = [176.919525856024,0.255696872807351]
        self.tranges = [[ 765579107.995, 765579821.995 ],
                        [ 766525332.995, 766526576.995 ],
                        [ 919754986.995, 919755916.995 ]]
        self.trange = tranges[0]
        self.parser = ga.setup_parser()
        self.args = self.parser.parse_args()

    def test_addhdr_default(self):
        self.assertFalse(args.addhdr)

    def test_annulus_default(self):
        self.assertIsNone(args.annulus)

    def test_annulus1_default(self):
        self.assertIsNone(args.annulus1)

    def test_annulus2_default(self):
        self.assertIsNone(args.annulus2)

    def test_band_default(self):
        self.assertEqual(args.band,'NUV')

    def test_calpath_default(self):
        self.assertEqual(args.calpath,'../cal/')

    def test_coadd_default(self):
        self.assertFalse(args.coadd)

    def test_csvfile_default(self):
        self.assertIsNone(args.csvfile)

    def test_dec_default(self):
        self.assertIsNone(args.dec)

    def test_detsize_default(self):
        self.assertAlmostEqual(args.detsize,1.25)

    def test_iocode_default(self):
        self.assertEqual(args.iocode,'wb')

    def test_maskdepth_default(self):
        self.assertAlmostEqual(args.maskdepth,20.0)

    def test_maskradius_default(self):
        self.assertAlmostEqual(args.maskradius=1.5)

    def test_maxgap_default(self):
        self.assertAlmostEqual(args.maxgap=1500.0

    def test_minexp_default(self):
        self.assertAlmostEqual(args.minexp,1.0))

    def test_overwrite_default(self):
        self.assertFalse(args.overwrite)

    def test_ra_default(self):
        self.assertIsNone(args.ra)

    def test_radius_default(self):
        self.assertIsNone(args.radius)

    def test_retries_default(self):
        self.assertEqual(args.retries,20)

    def test_skypos_default(self):
        self.assertIsNone(args.skypos)

    def test_stamp_default(self):
        self.assertIsNone(args.stamp)

    def test_stepsz_default(self):
        self.assertAlmostEqual(args.stepsz,0.0)

    def test_suggest_default(self):
        self.assertFalse(args.suggest)

    def test_tmax_default(self):
        self.assertAlmostEqual(args.tmax,1000000000000.0)

    def test_tmin_default(self):
        self.assertAlmostEqual(args.tmin,1.0)

    def test_trange_default(self):
        self.assertIsNone(args.trange)

    def test_verbose_default(self):
        self.assertEqual(args.verbose,0)

suite = unittest.TestLoader().loadTestsFromTestCase(TestgApertureArguments)
unittest.TextTestRunner(verbosity=2).run(suite)
