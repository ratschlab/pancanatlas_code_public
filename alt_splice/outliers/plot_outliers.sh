#!/bin/bash

set -e

### directory containing the results of find_outliers.py in tsv format
BASEDIR=
for event in $(grep -f ../annotation/census ${BASEDIR}/*.whitelisted.tsv | cut -f 2,5 | sort -u | tr '\t' '_')
do
   python sf_heatmap.py --event $event census 
done
