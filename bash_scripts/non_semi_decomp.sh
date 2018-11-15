#!/bin/bash

# Main script used for running reverse semi-decomposition algorithm.

cat ../dat_files/decomp_list | while read gnum
do 
    python3 ../main_sym.py $gnum > ../decomp/non_semidirect/decomp_${gnum}
done

