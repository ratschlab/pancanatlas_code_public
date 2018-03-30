#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/TCGA/PanCancer/rerun2018
datadir=/cluster/work/grlab/projects/TCGA/PanCancer/rerun2018_data
mem=80000
threads=8
pmem=$(($mem / $threads))

while IFS='' read -r line || [[ -n "$line" ]]
do
    aid=$(echo $line | cut -f 5 -d ' ')
    tid=$(echo $line | cut -f 3 -d ' ')
    pst=$(echo $line | cut -f 8 -d ' ')
    uuid=$(echo $line | cut -f 9 -d ' ')

    if [ "$pst" == "SINGLE" ]
    then    
        continue
    fi
    if [ "$uuid" == "NA" ]
    then
        continue
    fi
    if [ "$aid" == "aliquot_id" ]
    then
        continue
    fi
    donefile="${basedir}/${tid}.${aid}/pipeline.done"
    lockfile="${basedir}/${tid}.${aid}/running.lock"
    if [ -f $donefile ]
    then
        #echo $donefile exists
        continue
    elif [ ! -d "${datadir}/${tid}.${aid}" ]
    then
        echo ${tid}.${aid}: alignment not downloaded, yet.
    elif [ -z "$(find  ${datadir}/${tid}.${aid} -name \*.bam)" ]
    then
        echo ${tid}.${aid}: download incomplete
    elif [ -f $lockfile ]
    then
        echo $lockfile exists - job still running?
    else
        echo submitting ${tid} ${aid}
        echo "$(pwd)/process_sample_one.sh $tid $aid" | bsub -M ${mem} -W 12:00 -n ${threads} -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -J PCrr18 -o ${basedir}/${tid}.${aid}.cluster.log
    fi
done < /cluster/work/grlab/projects/TCGA/PanCancer/metadata/2018-02-03/full_metadata.overlap_RNA_WXS.RNA.GDC_ID.tsv
