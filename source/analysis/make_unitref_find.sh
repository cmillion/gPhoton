#!/bin/sh

# This script can be used to generate the StateFiles used as UnitTestReferenceFiles.  The output is not put directly into the UnitTestReferenceFiles folders to avoid accidental overwriting of the reference files.  If you need to replace the reference files for unit testing in the repo, run the appropriate case(s) and move them into position manually.

# Print starting time.
date

# This is used as a starting point for the various gtool_find tests.
python gtool_input.py --radeccoords 244.27246917 55.26919386 -i cr_dra --coordtype radec -o test_files/StateFiles/UnitTestReferenceFiles/gtool_find/ -n TestCase00a
# Get some exptimes populated that should be different from the final result of the unit tests below.
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00a_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00b_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00b_cr_dra.json -b BOTH --t0 806326407.995 --t1 868403745.995
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00a_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00c_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00c_cr_dra.json -b FUV --t0 806326407.995 --t1 868403745.995
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00a_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00d_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00d_cr_dra.json -b NUV --t0 806326407.995 --t1 868403745.995


# Test Case01: Band = FUV, Overwrite = Yes, Previous Exptime = No
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00a_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase01_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase01_cr_dra.json -b FUV

# Test Case02: Band = FUV, Overwrite = Yes, Previous Exptime = FUV-yes,NUV-yes
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00b_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase02_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase02_cr_dra.json -b FUV

# Test Case03: Band = FUV, Overwrite = Yes, Previous Exptime = FUV-yes,NUV-no
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00c_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase03_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase03_cr_dra.json -b FUV

# Test Case04: Band = FUV, Overwrite = Yes, Previous Exptime = FUV-no,NUV-yes
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00d_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase04_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase04_cr_dra.json -b FUV

# Test Case05: Band = FUV, Overwrite = No, Previous Exptime = No
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00a_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase05_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase05_cr_dra.json -b FUV --no_overwrite

# Test Case06: Band = FUV, Overwrite = No, Previous Exptime = FUV-yes,NUV-yes
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00b_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase06_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase06_cr_dra.json -b FUV --no_overwrite

# Test Case07: Band = FUV, Overwrite = No, Previous Exptime = FUV-yes,NUV-no
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00c_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase07_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase07_cr_dra.json -b FUV --no_overwrite

# Test Case08: Band = FUV, Overwrite = No, Previous Exptime = FUV-no,NUV-yes
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00d_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase08_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase08_cr_dra.json -b FUV --no_overwrite

# Test Case09: Band = NUV, Overwrite = Yes, Previous Exptime = No
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00a_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase09_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase09_cr_dra.json -b NUV

# Test Case10: Band = NUV, Overwrite = Yes, Previous Exptime = FUV-yes,NUV-yes
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00b_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase10_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase10_cr_dra.json -b NUV

# Test Case11: Band = NUV, Overwrite = Yes, Previous Exptime = FUV-yes,NUV-no
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00c_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase11_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase11_cr_dra.json -b NUV

# Test Case12: Band = NUV, Overwrite = Yes, Previous Exptime = FUV-no,NUV-yes
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00d_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase12_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase12_cr_dra.json -b NUV

# Test Case13: Band = NUV, Overwrite = No, Previous Exptime = No
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00a_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase13_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase13_cr_dra.json -b NUV --no_overwrite

# Test Case14: Band = NUV, Overwrite = No, Previous Exptime = FUV-yes,NUV-yes
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00b_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase14_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase14_cr_dra.json -b NUV --no_overwrite

# Test Case15: Band = NUV, Overwrite = No, Previous Exptime = FUV-yes,NUV-no
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00c_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase15_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase15_cr_dra.json -b NUV --no_overwrite

# Test Case16: Band = NUV, Overwrite = No, Previous Exptime = FUV-no,NUV-yes
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00d_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase16_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase16_cr_dra.json -b NUV --no_overwrite

# Test Case17: Band = BOTH, Overwrite = Yes, Previous Exptime = No
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00a_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase17_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase17_cr_dra.json -b BOTH

# Test Case18: Band = BOTH, Overwrite = Yes, Previous Exptime = FUV-yes,NUV-yes
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00b_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase18_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase18_cr_dra.json -b BOTH

# Test Case19: Band = BOTH, Overwrite = Yes, Previous Exptime = FUV-yes,NUV-no
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00c_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase19_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase19_cr_dra.json -b BOTH

# Test Case20: Band = BOTH, Overwrite = Yes, Previous Exptime = FUV-no,NUV-yes
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00d_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase20_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase20_cr_dra.json -b BOTH

# Test Case21: Band = BOTH, Overwrite = No, Previous Exptime = No
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00a_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase21_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase21_cr_dra.json -b BOTH --no_overwrite

# Test Case22: Band = BOTH, Overwrite = No, Previous Exptime = FUV-yes,NUV-yes
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00b_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase22_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase22_cr_dra.json -b BOTH --no_overwrite

# Test Case23: Band = BOTH, Overwrite = No, Previous Exptime = FUV-yes,NUV-no
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00c_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase23_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase23_cr_dra.json -b BOTH --no_overwrite

# Test Case24: Band = BOTH, Overwrite = No, Previous Exptime = FUV-no,NUV-yes
\cp -f test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase00d_cr_dra.json test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase24_cr_dra.json
python gtool_find.py test_files/StateFiles/UnitTestReferenceFiles/gtool_find/TestCase24_cr_dra.json -b BOTH --no_overwrite

# Print ending time.
date

