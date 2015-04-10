#!/bin/sh

if [ $# -ne 1 ]; then
    echo "Please choose a test to run by entering a number."
    exit 1
fi

if [ "$1" -eq 0 ]; then
    echo "Choose option {1,2,3,4,5,6,7,8}"
fi

if [ "$1" -eq 1 ]; then
    python gtool_input.py -f test_files/cr_dra.radec_noid.txt --coordtype radec -o test_files/StateFiles/
fi

if [ "$1" -eq 2 ]; then
    python gtool_input.py -f test_files/cr_dra.radec.txt --coordtype radec -o test_files/StateFiles/
fi

if [ "$1" -eq 3 ]; then
    python gtool_input.py -f test_files/cr_dra.galcoord_noid.txt --coordtype galactic -o test_files/StateFiles/
fi

if [ "$1" -eq 4 ]; then
    python gtool_input.py -f test_files/cr_dra.galcoord.txt --coordtype galactic -o test_files/StateFiles/
fi

if [ "$1" -eq 5 ]; then
    python gtool_input.py --radeccoords 244.27246917 55.26919386 -o test_files/StateFiles/
fi

if [ "$1" -eq 6 ]; then
    python gtool_input.py --radeccoords 244.27246917 55.26919386 -i cr_dra -o test_files/StateFiles/
fi

if [ "$1" -eq 7 ]; then
    python gtool_input.py --galcoords 84.90269369 43.70873321 -o test_files/StateFiles/
fi

if [ "$1" -eq 8 ]; then
    python gtool_input.py --galcoords 84.90269369 43.70873321 -i cr_dra -o test_files/StateFiles/
fi