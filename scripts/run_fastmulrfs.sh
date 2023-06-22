#!/bin/bash

input_trees=$1
temp_trees=$2
sptree=$3

preprocess="/gpfs/projects/bsc40/current/gmutti/software/phylogeny/fastmulrfs/python-tools/preprocess_multrees_v3.py"
fastmulrfs="/gpfs/projects/bsc40/current/gmutti/software/phylogeny/fastmulrfs/external/FastRFS/build/FastRFS"

python $preprocess \
    -i $input_trees \
    -o $temp_trees \
    --verbose

if [ ! -e $fastmulrfs ]; then
    echo "Need to get external dependencies!"
else
    $fastmulrfs \
        -i $temp_trees \
        -o $sptree 
fi
