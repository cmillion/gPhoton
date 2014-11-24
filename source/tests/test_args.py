import unittest
import gMap as gm
import gFind as gf
import gAperture as ga

"""Unit tests for the command line parameters and checking for the gPhoton
command line utilities [gPhoton.py, gAperture.py,gMap.py, gFind.py].

Most of the parameters propagate from gphoton_args to all three command line
utilities. We'll just test these through gAperture and then individually test
the few that don't propagate through gAperture.
"""
class TestArguments(unittest.TestCase):
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

    # BEGIN TESTING DEFAULTS
    def test_addhdr_default(self):
        """Check the default value of --addhdr"""
        self.assertFalse(self.args.addhdr)

    def test_annulus_default(self):
        """Check the default value of --annulus"""
        self.assertIsNone(self.args.annulus)

    def test_annulus1_default(self):
        """Check the default value of --annulus1"""
        self.assertIsNone(self.args.annulus1)

    def test_annulus2_default(self):
        """Check the default value of --annulus2"""
        self.assertIsNone(self.args.annulus2)

    def test_band_default(self):
        """Check the default value of --band"""
        self.assertEqual(self.args.band,'NUV')

    def test_calpath_default(self):
        """Check the default value of --calpath"""
        self.assertEqual(self.args.calpath,'../cal/')

    def test_coadd_default(self):
        """Check the default value of --coadd"""
        self.assertFalse(self.args.coadd)

    def test_csvfile_default(self):
        """Check the default value of --csvfile"""
        self.assertIsNone(self.args.csvfile)

    def test_dec_default(self):
        """Check the default value of --dec"""
        self.assertIsNone(self.args.dec)

    def test_detsize_default(self):
        """Check the default value of --detsize"""
        self.assertAlmostEqual(self.args.detsize,1.25)

    def test_iocode_default(self):
        """Check the default value of --iocode"""
        self.assertEqual(self.args.iocode,'wb')

    def test_maskdepth_default(self):
        """Check the default value of --bgmaskdepth"""
        self.assertAlmostEqual(self.args.maskdepth,20.0)

    def test_maskradius_default(self):
        """Check the default value of --bgmaskradius"""
        self.assertAlmostEqual(self.args.maskradius,1.5)

    def test_maxgap_default(self):
        """Check the default value of --maxgap"""
        self.assertAlmostEqual(self.args.maxgap,1500.0)

    def test_minexp_default(self):
        """Check the default value of --minexp"""
        self.assertAlmostEqual(self.args.minexp,1.0)

    def test_overwrite_default(self):
        """Check the default value of --overwrite"""
        self.assertFalse(self.args.overwrite)

    def test_ra_default(self):
        """Check the default value of --ra"""
        self.assertIsNone(self.args.ra)

    def test_radius_default(self):
        """Check the default value of --radius"""
        self.assertIsNone(self.args.radius)

    def test_retries_default(self):
        """Check the default value of --retries"""
        self.assertEqual(self.args.retries,20)

    def test_skypos_default(self):
        """Check the default value of --skypos"""
        self.assertIsNone(self.args.skypos)

    def test_stamp_default(self):
        """Check the default value of --stamp"""
        self.assertIsNone(self.args.stamp)

    def test_stepsz_default(self):
        """Check the default value of --stepsz"""
        self.assertAlmostEqual(self.args.stepsz,0.0)

    def test_suggest_default(self):
        """Check the default value of --suggest"""
        self.assertFalse(self.args.suggest)

    def test_tmax_default(self):
        """Check the default value of --tmax"""
        self.assertAlmostEqual(self.args.tmax,1000000000000.0)

    def test_tmin_default(self):
        """Check the default value of --tmin"""
        self.assertAlmostEqual(self.args.tmin,1.0)

    def test_trange_default(self):
        """Check the default value of --trange"""
        self.assertIsNone(self.args.trange)

    def test_verbose_default(self):
        """Check the default value of --verbose"""
        self.assertEqual(self.args.verbose,0)
    # END TESTING DEFAULTS

    def test_ra_dec(self):
        """Check that RA/Dec are recorded."""
        args = self.parser.parse_args(['--ra',str(self.skypos[0]),
            '--dec',str(self.skypos[1]),'-a',str(self.radius)])
        self.assertAlmostEqual(args.ra,self.skypos[0])
        self.assertAlmostEqual(args.dec,self.skypos[1])
        self.assertIsNone(args.skypos)

    def test_skypos_propagate(self):
        """Check that RA/Dec are propagated to skypos."""
        args = self.parser.parse_args(['--ra',str(self.skypos[0]),
            '--dec',str(self.skypos[1]),'-a',str(self.radius)])
        args = ga.check_args(args)
        self.assertAlmostEqual(args.skypos[0],self.skypos[0])
        self.assertAlmostEqual(args.skypos[1],self.skypos[1])

    def test_ra_check_toobig(self):
        """Check that RA is forced to have reasonable values"""
        args = self.parser.parse_args(['--ra',str(400),
            '--dec',str(self.skypos[1]),'-a',str(self.radius)])
        with self.assertRaises(SystemExit):
            ga.check_args(args)

    def test_dec_check_toobig(self):
        """Check that Dec is forced to have reasonable values"""
        args = self.parser.parse_args(['--ra',str(self.skypos[0]),
            '--dec',str(100),'-a',str(self.radius)])
        with self.assertRaises(SystemExit):
            ga.check_args(args)

    def test_ra_check_toosmall(self):
        """Check that RA is forced to have reasonable values"""
        args = self.parser.parse_args(['--ra',str(-20),
            '--dec',str(self.skypos[1]),'-a',str(self.radius)])
        with self.assertRaises(SystemExit):
            ga.check_args(args)

    def test_dec_check_toosmall(self):
        """Check that Dec is forced to have reasonable values"""
        args = self.parser.parse_args(['--ra',str(self.skypos[0]),
            '--dec',str(-100),'-a',str(self.radius)])
        with self.assertRaises(SystemExit):
            ga.check_args(args)

    def test_no_skypos1(self):
        """Check that either RA/Dec or skypos is specified"""
        args = self.parser.parse_args(['-a',str(self.radius)])
        with self.assertRaises(SystemExit):
            ga.check_args(args)

    def test_no_skypos2(self):
        """Check that either RA/Dec or skypos is specified"""
        args = self.parser.parse_args(['--ra',str(self.skypos[0])])
        with self.assertRaises(SystemExit):
            ga.check_args(args)

    def test_both_skypos(self):
        """Check that not _both_ RA/Dec and skypos are specified."""
        args = self.parser.parse_args(['--ra',str(self.skypos[0]),
            '--dec',str(self.skypos[1]),
            '--skypos',str([self.skypos[0]+1,self.skypos[1]])])
        with self.assertRaises(SystemExit):
            ga.check_args(args)

    def test_both_skypos_both(self):
        """Check that if both agreee, it's fine... weird edge case."""
        args = self.parser.parse_args(['--ra',str(self.skypos[0]),
            '--dec',str(self.skypos[1]),
            '--skypos',str([self.skypos[0],self.skypos[1]])])
        self.assertAlmostEqual(self.skypos[0],args.skypos[0])
        self.assertAlmostEqual(self.skypos[1],args.skypos[1])

    def test_no_radius(self):
        """Check that the aperture radius is mandatory."""
        args = self.parser.parse_args(['--skypos',str(self.skypos)])
        with self.assertRaises(SystemExit):
            ga.check_args(args)

    def test_annulus_propagate1(self):
        """Check that inner and outer get propagated to annulus."""
        args = self.parser.parse_args(['--skypos',str(self.skypos),
            '-a',str(self.radius),'-i',str(self.annulus[0]),
            '-o',str(self.annulus[1])])
        args = ga.check_args(args)
        self.assertAlmostEqual(self.annulus[0],args.annulus1)
        self.assertAlmostEqual(self.annulus[1],args.annulus2)
        self.assertAlmostEqual(self.annulus[0],args.annulus[0])
        self.assertAlmostEqual(self.annulus[1],args.annulus[1])

    def test_annulus_propagate2(self):
        """Check that annulus propagates to inner and outer."""
        args = self.parser.parse_args(['--skypos',str(self.skypos),
            '-a',str(self.radius),'--annulus',str(self.annulus)])
        args = ga.check_args(args)
        self.assertAlmostEqual(self.annulus[0],args.annulus1)
        self.assertAlmostEqual(self.annulus[1],args.annulus2)
        self.assertAlmostEqual(self.annulus[0],args.annulus[0])
        self.assertAlmostEqual(self.annulus[1],args.annulus[1])

    def test_tranges(self):
        self.assertTrue(True)

    def test_gaper_default(self):
        """Check the default value of --alt (gFind)"""
        self.assertFalse(gf.setup_parser().parse_args().gaper)

    def test_exponly_default(self):
        """Check the default value of --exponly (gFind)"""
        self.assertFalse(gf.setup_parser().parse_args().exponly)

    def test_quiet_default(self):
        """Check the default value of --quite (gFind)"""
        self.assertFalse(gf.setup_parser().parse_args().quiet)

    def test_angle_default(self):
        """Check the default value of --angle (gMap)"""
        self.assertIsNone(gm.setup_parser().parse_args().angle)

    def test_count_default(self):
        """Check the default value of --count (gMap)"""
        self.assertIsNone(gm.setup_parser().parse_args().cntfile)

    def test_raangle_default(self):
        """Check the default value of --raangle (gMap)"""
        self.assertIsNone(gm.setup_parser().parse_args().raangle)

    def test_decangle_default(self):
        """Check the default value of --decangle (gMap)"""
        self.assertIsNone(gm.setup_parser().parse_args().decangle)

    def test_intensity_default(self):
        """Check the default value of --intensity (gMap)"""
        self.assertIsNone(gm.setup_parser().parse_args().intfile)

    def test_memlight_default(self):
        """Check the default value of --memlight (gMap)"""
        self.assertAlmostEqual(100.,gm.setup_parser().parse_args().memlight)

    def test_response_default(self):
        """Check the default value of --response (gMap)"""
        self.assertIsNone(gm.setup_parser().parse_args().rrfile)

    def test_angle_propagate1(self):
        """Check that raangle and decangle propagate to skyrange (gMap)"""
        args = gm.check_args(gm.setup_parser().parse_args([
            '--skypos',str(self.skypos),'--raangle','0.1','--decangle','0.2']))
        self.assertAlmostEqual(0.1,args.skyrange[0])
        self.assertAlmostEqual(0.2,args.skyrange[1])

    def test_angle_propagate2(self):
        """Check that skyrange propagates to raangle and decangle. (gMap)"""
        args = gm.check_args(gm.setup_parser().parse_args(
            ['--skypos',str(self.skypos),'--skyrange','[0.1,0.2]']))
        self.assertAlmostEqual(0.1,args.raangle)
        self.assertAlmostEqual(0.2,args.decangle)

    def test_angle_both(self):
        """Check that both raangle and decangle must be specified. (gMap)"""
        with self.assertRaises(SystemExit):
            gm.check_args(gm.check_args(gm.setup_parser().parse_args(
                ['--skypos',str(self.skypos),'--raangle','0.1'])))

    def test_tmin_tmax_propagate_defaults(self):
        """If trange is not given, sub in tmin/tmax."""
        args = self.parser.parse_args(['--skypos',str(self.skypos),
            '-a',str(self.radius)])
        args = ga.check_args(args)
        self.assertAlmostEqual(args.tmin,args.trange[0][0])
        self.assertAlmostEqual(args.tmax,args.trange[0][1])

suite = unittest.TestLoader().loadTestsFromTestCase(TestArguments)
unittest.TextTestRunner(verbosity=2).run(suite)
