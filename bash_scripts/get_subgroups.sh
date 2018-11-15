#!/bin/bash

# Script used to download Bieberbach and symmorphic space subgroups for each space group.

cat ../dat_files/groups_to_decompose | while read gnum
do
    python3 ../get_subgroups.py $gnum 
done

