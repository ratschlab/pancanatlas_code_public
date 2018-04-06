#!/bin/bash

set -e

# path to annotation in GTF format
anno=
# path to the spladder executable
as_script=
# text file containing the absolute paths to the alignment files in BAM format
bam_files= 
# directory for the output files
outdir= 

### run spladder 
python ${as_script} -b ${bam_files} -o $outdir -a ${anno} -v y -c 2 -M merge_graphs -T y -V y -n 50 -P y -p y --sparse_bam y -D 500 -t exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons
