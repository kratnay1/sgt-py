#!/bin/bash

# Script previously used to create table entries for reverse semi-decompositions.

dir=../tex

ls $dir/table_entries/non_semidirect/entry_?? | while read filename
do 
    cat $filename >> $dir/table_sym.tex
done

ls $dir/table_entries/non_semidirect/entry_??? | while read filename
do 
    cat $filename >> $dir/table_sym.tex
done

