#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <type>"
    echo "  type is one of alt or non-alt"
    exit 1
else
    typ="$1"
fi

basedir=/cluster/work/grlab/projects/TCGA/PanCancer/rerun2018
outdir=${basedir}_hdf5
mkdir -p $outdir
if [ "$typ" == "alt" ]
then
    outfile=${outdir}/expression_counts.hdf5
    python collect_counts_incrementally.py $outfile "${basedir}/*/*.aligned.ex_cnt.tab" -v
elif [ "$typ" == "non-alt" ]
then
    outfile=${outdir}/expression_counts.non_alt.hdf5
    python collect_counts_incrementally.py $outfile "${basedir}/*/*.aligned.ex_cnt.non_alt.tab" -v
else
    echo "Unknown type: $typ"
    exit 1
fi
