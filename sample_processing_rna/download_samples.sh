#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/TCGA/PanCancer/rerun2018_data
resultdir=/cluster/work/grlab/projects/TCGA/PanCancer/rerun2018
gdc=<PATH TO GDC CLIENT>
keys=<FOLDER with KEYS>

mkdir -p $basedir
cnt=0
cnt_done=0
while IFS='' read -r line || [[ -n "$line" ]]
do
    aid=$(echo $line | cut -f 5 -d ' ')
    tid=$(echo $line | cut -f 3 -d ' ')
    pair=$(echo $line | cut -f 8 -d ' ')
    uuid=$(echo $line | cut -f 9 -d ' ')
    if [ "$aid" == "aliquot_id" ]
    then 
        continue
    fi
    if [ "$uuid" == "NA" ]
    then
        continue
    fi
    if [ "$pair" != "PAIRED" ]
    then
        continue
    fi

    ### clean up on the way
    outdir=${basedir}/${tid}.${aid}
    if [ -f "${resultdir}/${tid}.${aid}/pipeline.done" ]
    then
        cnt_done=$(($cnt_done + 1))
        continue
    fi

    mkdir -p $outdir
    #for uuid in $(python get_uuid.py $aid)
    #do  

    if [ ! -d ${outdir}/${uuid} ]
    then
        echo "Downloading to $outdir $uuid"
        $gdc download -t ${keys}/<REDACTED> -n 2 -d $outdir $uuid &
        cnt=$(($cnt + 1))
    else
        echo " ${outdir}/${uuid} already exists"
    fi
    if [ "$cnt" == "8" ]
    then
        wait
        cnt=0
    fi

    #done
done < full_metadata.overlap_RNA_WXS.RNA.GDC_ID.tsv
echo $cnt_done samples were already finished
