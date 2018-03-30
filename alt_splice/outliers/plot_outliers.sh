#!/bin/bash

set -e

BASEDIR="/cluster/work/grlab/projects/TCGA/PanCancer/rerun2018_alt_splice/outliers/"
for event in $(grep -f ${BASEDIR}/census ${BASEDIR}/*.whitelisted.tsv | cut -f 2,5 | sort -u | tr '\t' '_')
do
   python sf_heatmap.py --event $event census 
done
