__version__ = '1.0'

"""
.. module:: test_gtool_input

   :synopsis: Test module for gtool_input.py

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import os
import gtool_input
import unittest

#--------------------
class TestgToolInput(unittest.TestCase):
    """ Note: Coordinates are passed as strings to replicate 
    command-line usage. """

    reference_file_path = ("test_files/StateFiles/"
                           "UnitTestReferenceFiles/")
    def testCase01(self):
        output_files = gtool_input.input_targets(ifile = 
                                  "test_files/cr_dra.radec_noid.txt",
                                  coordtype = "radec", 
                                  output_dir = "test_files/StateFiles/",
                                  output_prefix = "TestCase01")
        for f in output_files:
            with open(f, 'rb') as nf:
                new_content = nf.readlines()
            old_file = (self.reference_file_path + "Case01/" + 
                      os.path.basename(f))
            if os.path.isfile(old_file):
                with open(old_file, 'rb') as of:
                    old_content = of.readlines()
            else:
                self.fail(msg="Reference file not found.  Looking for"
                          " "+old_file)
            self.assertEqual(old_content, new_content)
            os.remove(f)

    def testCase02(self):
        output_files = gtool_input.input_targets(ifile = 
                                  "test_files/cr_dra.radec.txt",
                                  coordtype = "radec", 
                                  output_dir = "test_files/StateFiles/",
                                  output_prefix = "TestCase02")
        for f in output_files:
            with open(f, 'rb') as nf:
                new_content = nf.readlines()
            old_file = (self.reference_file_path + "Case02/" + 
                      os.path.basename(f))
            if os.path.isfile(old_file):
                with open(old_file, 'rb') as of:
                    old_content = of.readlines()
            else:
                self.fail(msg="Reference file not found.  Looking for"
                          " "+old_file)
            self.assertEqual(old_content, new_content)
            os.remove(f)

    def testCase03(self):
        output_files = gtool_input.input_targets(ifile = 
                                  "test_files/cr_dra.galcoord_noid.txt",
                                  coordtype = "galactic", 
                                  output_dir = "test_files/StateFiles/",
                                  output_prefix = "TestCase03")
        for f in output_files:
            with open(f, 'rb') as nf:
                new_content = nf.readlines()
            old_file = (self.reference_file_path + "Case03/" + 
                      os.path.basename(f))
            if os.path.isfile(old_file):
                with open(old_file, 'rb') as of:
                    old_content = of.readlines()
            else:
                self.fail(msg="Reference file not found.  Looking for"
                          " "+old_file)
            self.assertEqual(old_content, new_content)
            os.remove(f)

    def testCase04(self):
        output_files = gtool_input.input_targets(ifile = 
                                  "test_files/cr_dra.galcoord.txt",
                                  coordtype = "galactic", 
                                  output_dir = "test_files/StateFiles/",
                                  output_prefix = "TestCase04")
        for f in output_files:
            with open(f, 'rb') as nf:
                new_content = nf.readlines()
            old_file = (self.reference_file_path + "Case04/" + 
                      os.path.basename(f))
            if os.path.isfile(old_file):
                with open(old_file, 'rb') as of:
                    old_content = of.readlines()
            else:
                self.fail(msg="Reference file not found.  Looking for"
                          " "+old_file)
            self.assertEqual(old_content, new_content)
            os.remove(f)

    def testCase05(self):
        output_files = gtool_input.input_targets(radeccoords = 
                                                 ['244.27246917',
                                                 '55.26919386'],
                                  output_dir = "test_files/StateFiles/",
                                  output_prefix = "TestCase05")
        for f in output_files:
            with open(f, 'rb') as nf:
                new_content = nf.readlines()
            old_file = (self.reference_file_path + "Case05/" + 
                      os.path.basename(f))
            if os.path.isfile(old_file):
                with open(old_file, 'rb') as of:
                    old_content = of.readlines()
            else:
                self.fail(msg="Reference file not found.  Looking for"
                          " "+old_file)
            self.assertEqual(old_content, new_content)
            os.remove(f)

    def testCase06(self):
        output_files = gtool_input.input_targets(radeccoords = 
                                                 ['244.27246917',
                                                 '55.26919386'],
                                  target_id = "cr_dra",
                                  output_dir = "test_files/StateFiles/",
                                  output_prefix = "TestCase06")
        for f in output_files:
            with open(f, 'rb') as nf:
                new_content = nf.readlines()
            old_file = (self.reference_file_path + "Case06/" + 
                      os.path.basename(f))
            if os.path.isfile(old_file):
                with open(old_file, 'rb') as of:
                    old_content = of.readlines()
            else:
                self.fail(msg="Reference file not found.  Looking for"
                          " "+old_file)
            self.assertEqual(old_content, new_content)
            os.remove(f)

    def testCase07(self):
        output_files = gtool_input.input_targets(galcoords = 
                                                 ['84.90269369',
                                                 '43.70873321'],
                                  output_dir = "test_files/StateFiles/",
                                  output_prefix = "TestCase07")
        for f in output_files:
            with open(f, 'rb') as nf:
                new_content = nf.readlines()
            old_file = (self.reference_file_path + "Case07/" + 
                      os.path.basename(f))
            if os.path.isfile(old_file):
                with open(old_file, 'rb') as of:
                    old_content = of.readlines()
            else:
                self.fail(msg="Reference file not found.  Looking for"
                          " "+old_file)
            self.assertEqual(old_content, new_content)
            os.remove(f)

    def testCase08(self):
        output_files = gtool_input.input_targets(galcoords = 
                                                 ['84.90269369',
                                                 '43.70873321'],
                                  target_id = "cr_dra",
                                  output_dir = "test_files/StateFiles/",
                                  output_prefix = "TestCase08")
        for f in output_files:
            with open(f, 'rb') as nf:
                new_content = nf.readlines()
            old_file = (self.reference_file_path + "Case08/" + 
                      os.path.basename(f))
            if os.path.isfile(old_file):
                with open(old_file, 'rb') as of:
                    old_content = of.readlines()
            else:
                self.fail(msg="Reference file not found.  Looking for"
                          " "+old_file)
            self.assertEqual(old_content, new_content)
            os.remove(f)
    
#--------------------

if __name__ == "__main__":
    unittest.main()
