#!/bin/bash

set -e

### directory containing the results of the SplAdder run
basedir=
CONF=2

for file in ${basedir}/*C${CONF}.counts.hdf5 
do
    if [ -f "${file%hdf5}psi_filt_per_ct.pickle" ]
    then
        continue   
    fi

    echo "handling $file"
    if [ "$1" == "local" ]
    then
        python $(pwd)/gen_psi_filter_counts_per_ct.py $file
    else
        echo "python $(pwd)/gen_psi_filter_counts_per_ct.py $file" | bsub -M 10000 -J filt -n 1 -W 12:00 -R "rusage[mem=10000]" -R "span[hosts=1]" -o /dev/null
    fi
done
