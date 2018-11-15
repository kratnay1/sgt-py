#!/bin/bash

# Script previously used to create table entries for semi-decompositions.

dir=../tex

ls $dir/table_entries/semidirect/entry_?? | while read filename
do 
    cat $filename >> $dir/table.tex
done

ls $dir/table_entries/semidirect/entry_??? | while read filename
do 
    cat $filename >> $dir/table.tex
done

