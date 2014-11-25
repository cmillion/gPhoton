import unittest
import gMap as gm
import gFind as gf
import gAperture as ga

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
        self.parser = ga.setup_parser()
        self.args = self.parser.parse_args()

    def test_basic_query(self):
        """Test a simple query."""
        out = gf.gFind(band=self.bands[0],skypos=self.skypos)
        self.assertAlmostEqual(out[self.bands[0]]['expt'],2885.0)
        self.assertEqual(len(out[self.bands[0]]['t0']),4)

    def test_maxgap_query(self):
        """Test that the maxgap keyword combines nearby time ranges."""
        out = gf.gFind(band=self.bands[0],skypos=self.skypos,maxgap=100)
        self.assertAlmostEqual(out[self.bands[0]]['expt'],2888.0)
        self.assertEqual(len(out[self.bands[0]]['t0']),3)

    def test_minexp_query(self):
        """Test that the minexp keyword excludes too-short time ranges."""
        out = gf.gFind(band=self.bands[0],skypos=self.skypos,minexp=1000)
        self.assertAlmostEqual(out[self.bands[0]]['expt'],1244.0)
        self.assertEqual(len(out[self.bands[0]]['t0']),1)

    def test_detsize_query(self):
        """Test that the detsize keyword excludes exposure times."""
        out = gf.gFind(band=self.bands[0],skypos=self.skypos,detsize=0.5)
        self.assertAlmostEqual(out[self.bands[0]]['expt'],1506.0)
        self.assertEqual(len(out[self.bands[0]]['t0']),25)

    def test_both_query(self):
        """Test that band='both' does, in fact, search on both bands."""
        fuvout = gf.gFind(band='FUV',skypos=self.skypos)
        nuvout = gf.gFind(band='NUV',skypos=self.skypos)
        bothout = gf.gFind(band='both',skypos=self.skypos)
        self.assertAlmostEqual(fuvout['FUV']['expt'],bothout['FUV']['expt'])
        self.assertAlmostEqual(nuvout['NUV']['expt'],bothout['NUV']['expt'])
        self.assertEqual(len(fuvout['FUV']['t0']),len(bothout['FUV']['t0']))
        self.assertEqual(len(nuvout['NUV']['t0']),len(bothout['NUV']['t0']))

suite = unittest.TestLoader().loadTestsFromTestCase(TestRegression)
unittest.TextTestRunner(verbosity=2).run(suite)
