import unittest
import gPhoton.gMap as gm
import gPhoton.gFind as gf
import gPhoton.gAperture as ga

from gPhoton.gFind import setup_parser as gf_setup_parser

"""Regression tests for the command line utilities.
Note that these are end-to-end tests of the system that may fail due to
either changes in the source code or the database.
"""
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
        self.parser = gf_setup_parser()
        self.args = self.parser.parse_args()

    def test_basic_query(self):
        """Test a simple query."""
        out = gf(band=self.bands[0],skypos=self.skypos,quiet=True)
        self.assertAlmostEqual(out[self.bands[0]]['expt'],3516.0)
        self.assertEqual(len(out[self.bands[0]]['t0']),6)

    def test_maxgap_query(self):
        """Test that the maxgap keyword combines nearby time ranges."""
        out = gf(band=self.bands[0],skypos=self.skypos,maxgap=6000,
                       quiet=True)
        self.assertAlmostEqual(out[self.bands[0]]['expt'],9031.0)
        self.assertEqual(len(out[self.bands[0]]['t0']),5)

    def test_minexp_query(self):
        """Test that the minexp keyword excludes too-short time ranges."""
        out = gf(band=self.bands[0],skypos=self.skypos,minexp=1000,
                       quiet=True)
        self.assertAlmostEqual(out[self.bands[0]]['expt'],1244.0)
        self.assertEqual(len(out[self.bands[0]]['t0']),1)

    def test_detsize_query(self):
        """Test that the detsize keyword excludes exposure times."""
        out = gf(band=self.bands[0],skypos=self.skypos,detsize=0.5,
                       quiet=True)
        self.assertAlmostEqual(out[self.bands[0]]['expt'],2781.0)
        self.assertEqual(len(out[self.bands[0]]['t0']),5)

    def test_both_query(self):
        """Test that band='both' does, in fact, search on both bands."""
        fuvout = gf(band='FUV',skypos=self.skypos,quiet=True)
        nuvout = gf(band='NUV',skypos=self.skypos,quiet=True)
        bothout = gf(band='both',skypos=self.skypos,quiet=True)
        self.assertAlmostEqual(fuvout['FUV']['expt'],bothout['FUV']['expt'])
        self.assertAlmostEqual(nuvout['NUV']['expt'],bothout['NUV']['expt'])
        self.assertEqual(len(fuvout['FUV']['t0']),len(bothout['FUV']['t0']))
        self.assertEqual(len(nuvout['NUV']['t0']),len(bothout['NUV']['t0']))

suite = unittest.TestLoader().loadTestsFromTestCase(TestRegression)
unittest.TextTestRunner(verbosity=2).run(suite)
