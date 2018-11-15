#!/bin/bash

# Script needed to create table entries for groups where no decompositions were found.

cat ../dat_files/no_decomps.dat | while read gnum
do 
    python3 ../results.py $gnum
done

