#!/bin/sh

# This script can be used to generate the StateFiles used as UnitTestReferenceFiles.  The output is not put directly into the UnitTestReferenceFiles folders to avoid accidental overwriting of the reference files.  If you need to replace the reference files for unit testing in the repo, run the appropriate case(s) and move them into position manually.

# Test Case01: FromFile=YES, Coord=RADEC, ID=NO
python gtool_input.py -f test_files/cr_dra.radec_noid.txt --coordtype radec -o test_files/StateFiles/UnitTestReferenceFiles/gtool_input/ -n TestCase01

# Test Case02: FromFile=YES, Coord=RADEC, ID=YES
python gtool_input.py -f test_files/cr_dra.radec.txt --coordtype radec -o test_files/StateFiles/UnitTestReferenceFiles/gtool_input/ -n TestCase02

# Test Case03: FromFile=YES, Coord=GALACTIC, ID=NO
python gtool_input.py -f test_files/cr_dra.galcoord_noid.txt --coordtype galactic -o test_files/StateFiles/UnitTestReferenceFiles/gtool_input/ -n TestCase03

# Test Case04: FromFile=YES, Coord=GALACTIC, ID=YES
python gtool_input.py -f test_files/cr_dra.galcoord.txt --coordtype galactic -o test_files/StateFiles/UnitTestReferenceFiles/gtool_input/ -n TestCase04

# Test Case05: FromFile=NO, Coord=RADEC, ID=NO
python gtool_input.py --radeccoords 244.27246917 55.26919386 -o test_files/StateFiles/UnitTestReferenceFiles/gtool_input/ -n TestCase05

# Test Case06: FromFile=NO, Coord=RADEC, ID=YES
python gtool_input.py --radeccoords 244.27246917 55.26919386 -i cr_dra -o test_files/StateFiles/UnitTestReferenceFiles/gtool_input/ -n TestCase06

# Test Case07: FromFile=NO, Coord=GALACTIC, ID=NO
python gtool_input.py --galcoords 84.90269369 43.70873321 -o test_files/StateFiles/UnitTestReferenceFiles/gtool_input/ -n TestCase07

# Test Case08: FromFile=NO, Coord=GALACTIC, ID=YES
python gtool_input.py --galcoords 84.90269369 43.70873321 -i cr_dra -o test_files/StateFiles/UnitTestReferenceFiles/gtool_input/ -n TestCase08
