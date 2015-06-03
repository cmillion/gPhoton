__version__ = '1.0'

"""
.. module:: test_gtool_find

   :synopsis: Test module for gtool_find.py

.. moduleauthor:: Scott W. Fleming <fleming@stsci.edu>
"""

import os
import gtool_find
import unittest

#--------------------
class TestgToolFind(unittest.TestCase):

    reference_file_path = ("test_files/StateFiles/"
                           "UnitTestReferenceFiles/gtool_find/")
    def testCase01(self):
        pass

#--------------------

if __name__ == "__main__":
    unittest.main()
