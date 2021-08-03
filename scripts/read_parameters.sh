#!/bin/bash

for filename in ./*.txt; do
    echo $filename > params.txt
    echo `head -n 1 $filename` > params.txt
    echo > params.txt
done
