#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/TCGA/PanCancer/rerun2018_alt_splice_sub
CONF=2

for et in exon_skip intron_retention alt_3prime alt_5prime
do
    echo "handling $et"
    if [ "$1" == "local" ]
    then
        python $(pwd)/compare_events_anno.gtex_subset.py $et
    else
        echo "python $(pwd)/compare_events_anno.gtex_subset.py $et" | bsub -M 4000 -J filt -n 1 -W 4:00 -R "rusage[mem=4000]" -R "span[hosts=1]" -o /dev/null
    fi
done
