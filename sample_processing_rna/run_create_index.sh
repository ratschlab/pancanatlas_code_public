#!/bin/bash

set -e

star=/cluster/work/grlab/share/modules/packages/star/2.5.3a/bin/STAR
outdir=/cluster/work/grlab/projects/TCGA/PanCancer/genome/hg19_hs37d5_G19.overhang100_STAR_rr18
genomeFasta=/cluster/work/grlab/projects/TCGA/PanCancer/genome/hg19_hs37d5/genome.fa
annotation=/cluster/work/grlab/projects/TCGA/PanCancer/annotation/gencode.v19.annotation.hs37d5_chr.gtf
$star --runMode genomeGenerate --genomeDir $outdir --genomeFastaFiles $genomeFasta --sjdbOverhang 100 --sjdbGTFfile $annotation --runThreadN 4 > ${outdir}/run_star_index.log
