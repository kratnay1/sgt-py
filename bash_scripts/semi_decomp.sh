#!/bin/bash

# Main script used for running semi-decomposition algorithm.

cat ../dat_files/decomp_list | while read gnum
do 
    python3 ../main.py $gnum > ../decomp/semidirect/decomp_${gnum}
done

