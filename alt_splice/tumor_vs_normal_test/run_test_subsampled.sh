#!/bin/bash

set -e 

### absolute path to spladder_test.py
spladder=

### directory containing the results of the individual tests
basedir=

threads=4
CONF=2

studies="BLCA BRCA COAD HNSC KIRC KIRP LIHC LUAD LUSC PRAD READ STAD THCA UCEC" 
for run in $(seq 1 10)
do
    metadir=${basedir}/sample_lists/run_${run}
    for study in $studies
    do
        for event_type in alt_3prime alt_5prime exon_skip intron_retention mutex_exons
        do
            echo running $event_type on study $study
            outdir=${basedir}/testing_${study}-T_vs_${study}-N_R${run}
            if [ -f ${outdir}/test_results_C${CONF}_${event_type}.gene_unique.tsv ]
            then
                echo "already done"
            else
                mkdir -p $outdir
                echo "python $spladder -c $CONF -o $basedir -a ${metadir}/${study}_tumor.txt -b ${metadir}/${study}_normal.txt --out-tag R${run} --labelA ${study}-T --labelB ${study}-N -v y --event_types ${event_type} -v y --parallel $threads -D y -0 0.1 -V y" | bsub -M 200000 -J difftest -n $threads -R "rusage[mem=50000]" -R "span[hosts=1]" -o ${outdir}/cluster_${event_type}.log -W 24:00
            fi
        done
    done
done
