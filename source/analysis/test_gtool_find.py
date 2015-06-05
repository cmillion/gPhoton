__version__ = '1.0'

"""
.. module:: test_gtool_find

   :synopsis: Test module for gtool_find.py

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import os
import shutil
import gtool_find
import unittest
if __package__ is None:
    import sys, os
    gtfind_dir = os.path.dirname(os.path.abspath(__file__))
    gtfind_pardir = os.path.dirname(gtfind_dir)
    sys.path.insert(1, gtfind_pardir)
    from gFind import (setup_parser as gf_setup_parser)


#--------------------
class TestgToolFind(unittest.TestCase):

    """ Create default argument options. """
    parser = gtool_find.setup_args()
    parser = gf_setup_parser(parser=parser)
    args = vars(parser.parse_args())

    """ Define the base path to the unit test files. """
    reference_file_path = ("test_files/StateFiles/"
                           "UnitTestReferenceFiles/gtool_find/")
    def testCase01(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00a_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase01_cr_dra.json')
        ref_file = (self.reference_file_path + "Case01/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "FUV"
        self.args["no_overwrite"] = False
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase02(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00b_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase02_cr_dra.json')
        ref_file = (self.reference_file_path + "Case02/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "FUV"
        self.args["no_overwrite"] = False
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase03(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00c_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase03_cr_dra.json')
        ref_file = (self.reference_file_path + "Case03/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "FUV"
        self.args["no_overwrite"] = False
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase04(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00d_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase04_cr_dra.json')
        ref_file = (self.reference_file_path + "Case04/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "FUV"
        self.args["no_overwrite"] = False
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase05(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00a_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase05_cr_dra.json')
        ref_file = (self.reference_file_path + "Case05/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "FUV"
        self.args["no_overwrite"] = True
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase06(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00b_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase06_cr_dra.json')
        ref_file = (self.reference_file_path + "Case06/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "FUV"
        self.args["no_overwrite"] = True
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase07(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00c_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase07_cr_dra.json')
        ref_file = (self.reference_file_path + "Case07/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "FUV"
        self.args["no_overwrite"] = True
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase08(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00d_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase08_cr_dra.json')
        ref_file = (self.reference_file_path + "Case08/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "FUV"
        self.args["no_overwrite"] = True
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase09(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00a_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase09_cr_dra.json')
        ref_file = (self.reference_file_path + "Case09/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "NUV"
        self.args["no_overwrite"] = False
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase10(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00b_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase10_cr_dra.json')
        ref_file = (self.reference_file_path + "Case10/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "NUV"
        self.args["no_overwrite"] = False
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase11(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00c_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase11_cr_dra.json')
        ref_file = (self.reference_file_path + "Case11/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "NUV"
        self.args["no_overwrite"] = False
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase12(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00d_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase12_cr_dra.json')
        ref_file = (self.reference_file_path + "Case12/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "NUV"
        self.args["no_overwrite"] = False
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase13(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00a_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase13_cr_dra.json')
        ref_file = (self.reference_file_path + "Case13/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "NUV"
        self.args["no_overwrite"] = True
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase14(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00b_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase14_cr_dra.json')
        ref_file = (self.reference_file_path + "Case14/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "NUV"
        self.args["no_overwrite"] = True
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase15(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00c_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase15_cr_dra.json')
        ref_file = (self.reference_file_path + "Case15/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "NUV"
        self.args["no_overwrite"] = True
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase16(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00d_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase16_cr_dra.json')
        ref_file = (self.reference_file_path + "Case16/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "NUV"
        self.args["no_overwrite"] = True
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase17(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00a_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase17_cr_dra.json')
        ref_file = (self.reference_file_path + "Case17/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "BOTH"
        self.args["no_overwrite"] = False
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase18(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00b_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase18_cr_dra.json')
        ref_file = (self.reference_file_path + "Case18/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "BOTH"
        self.args["no_overwrite"] = False
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase19(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00c_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase19_cr_dra.json')
        ref_file = (self.reference_file_path + "Case19/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "BOTH"
        self.args["no_overwrite"] = False
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase20(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00d_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase20_cr_dra.json')
        ref_file = (self.reference_file_path + "Case20/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "BOTH"
        self.args["no_overwrite"] = False
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase21(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00a_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase21_cr_dra.json')
        ref_file = (self.reference_file_path + "Case21/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "BOTH"
        self.args["no_overwrite"] = True
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase22(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00b_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase22_cr_dra.json')
        ref_file = (self.reference_file_path + "Case22/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "BOTH"
        self.args["no_overwrite"] = True
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase23(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00c_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase23_cr_dra.json')
        ref_file = (self.reference_file_path + "Case23/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "BOTH"
        self.args["no_overwrite"] = True
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

    def testCase24(self):
        src_file = (self.reference_file_path + 
                    "Case00/TestCase00d_cr_dra.json")
        dest_file = (self.reference_file_path + 
                     'TestCase24_cr_dra.json')
        ref_file = (self.reference_file_path + "Case24/" + 
                    os.path.basename(dest_file))
        shutil.copyfile(src_file, dest_file)
        self.args["band"] = "BOTH"
        self.args["no_overwrite"] = True
        gtool_find.gtool_find(self.args, ifile=dest_file)
        with open(dest_file, 'rb') as nf:
            new_content = nf.readlines()
        if os.path.isfile(ref_file):
            with open(ref_file, 'rb') as of:
                old_content = of.readlines()
        else:
            self.fail(msg="Reference file not found.  Looking for " 
                      + ref_file)
        self.assertEqual(old_content, new_content)
        os.remove(dest_file)

#--------------------

if __name__ == "__main__":
    unittest.main()
