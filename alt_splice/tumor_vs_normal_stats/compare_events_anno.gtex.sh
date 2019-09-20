#!/bin/bash

set -e

### directory containing the results of the SplAdder run 
basedir=
CONF=2

for et in exon_skip intron_retention alt_3prime alt_5prime mutex_exons
do
    echo "handling $et"
    if [ "$1" == "local" ]
    then
        python $(pwd)/compare_events_anno.gtex.py $et
    else
        echo "python $(pwd)/compare_events_anno.gtex.py $et" | bsub -M 4000 -J filt -n 1 -W 4:00 -R "rusage[mem=4000]" -R "span[hosts=1]" -o /dev/null
    fi
done
