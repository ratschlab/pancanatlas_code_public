#!/bin/bash

set -e

module load jdk/8u92 

indir=$1
outdir=$2

script=/cluster/home/akahles/software/fastqc/FastQC_v0.11.6/fastqc

failed=N
any=N
for ffile in $(find $indir -name \*.fq)
do
    any=T
    $script -t 4 -o $outdir $ffile > ${outdir}/$(basename $ffile).fastqc.log 2>&1 || failed=T
done

if [ "$failed" == N -a "$any" == T ]
then
    touch "${outdir}/fastqc.done"
fi
