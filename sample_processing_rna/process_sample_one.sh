#!/bin/bash

set -e

### This script is the general processing wrapper for a single
### sample to be locall processed on a cluster node. The steps
### covered by this wrapper are:
### - fastqc
### - alignment & stats
### - expression counting
### - alternative splicing analysis (stage 1)


if [ -z "$2" ]
then
    echo "Usage: $0 <tcga_id> <aliquot_id>"
    exit 1
fi
tid=$1
aid=$2


### PROLOGUE ###
################

### software and fixed data patps
fastqc_script=${HOME}/git/projects/2013/PanCancerTCGA/rerun2017/sample_processing_rna/run_fastqc_one.sh
stats_script=${HOME}/git/projects/2013/PanCancerTCGA/rerun2015/sample_processing_rna/collect_quick_align_stats.py
count_script=${HOME}/git/projects/2013/PanCancerTCGA/count_expression/count_expression.py
as_script=${HOME}/git/software/spladder_pancan_rerun2018/python/spladder.py
bam2fastq=/${HOME}/git/projects/2013/PanCancerTCGA/rerun2017/sample_processing_rna/bam2fastq.py

samtools=/cluster/work/grlab/share/modules/packages/samtools/2.6/bin/samtools
stardir=/cluster/work/grlab/share/modules/packages/star/2.5.3a/bin
genome=/cluster/work/grlab/projects/TCGA/PanCancer/genome/hg19_hs37d5_G19.overhang100_STAR_rr18
genomeFasta=/cluster/work/grlab/projects/TCGA/PanCancer/genome/hg19_hs37d5/genome.fa
annotation=/cluster/work/grlab/projects/TCGA/PanCancer/annotation/gencode.v19.annotation.hs37d5_chr.gtf
annotation_spladder=/cluster/work/grlab/projects/TCGA/PanCancer/annotation/gencode.v19.annotation.hs37d5_chr.coding.spladder_rr18.gtf
data_dir=/cluster/work/grlab/projects/TCGA/PanCancer/rerun2018_data

threads=6
result_dir=/cluster/work/grlab/projects/TCGA/PanCancer/rerun2018/${tid}.${aid}
mkdir -p $result_dir
touch ${result_dir}/running.lock

export PATH=$stardir:$PATH

### create log dir
logdir=${result_dir}/logs
mkdir -p $logdir

### establish local work dir
workdir=$(mktemp -d -p $TMPDIR)
mkdir -p ${workdir}

echo "$(date) -- Starting analysis for $tid $aid"
echo ""
echo "$(date) -- Running with following configuration:"
echo "hostname: $(hostname)"
echo "result_dir: $result_dir"
echo "work_dir: $workdir"
echo "threads: $threads"
echo ""

### PREPROCESS BAM ###
######################
echo "$(date) -- Starting pre-processing of $tid $aid"
mkdir -p ${workdir}/orig
if [ -z "$(find ${data_dir}/${tid}.${aid} -name \*.bam)" ]
then
    echo "No alignment file for $tid $aid could be found"
    exit 1
fi
### symlink bam files
for file in $(find ${data_dir}/${tid}.${aid} -name \*.bam)
do
    filebase=$(basename $file)
    ln -s $file ${workdir}/orig/${filebase}
    ### sort bam files
    $samtools sort -n ${workdir}/orig/${filebase} -O BAM -T $TMPDIR --threads 3 > ${workdir}/orig/${filebase%bam}sorted.bam
    ### translate to fastq
    python $bam2fastq ${workdir}/orig/${filebase%bam}sorted.bam 50
done
echo "$(date) -- Pre-processing done"

### QC ###
##########
echo "$(date) -- Running QC in background"
mkdir -p ${workdir}/fastqc
echo "$fastqc_script ${workdir}/orig/ ${workdir}/fastqc "
$fastqc_script ${workdir}/orig/ ${workdir}/fastqc &


### ALIGNMENT ###
#################
echo "$(date) -- Starting alignment of $tid $aid"
mkdir -p ${workdir}/align
outfile_aln="${workdir}/align/${tid}.${aid}.aligned.bam"
outfile_jct="${workdir}/align/${tid}.${aid}.aligned.junctions"
donefile_aln="${workdir}/align/${tid}.${aid}.done"

### run STAR wrapper
cd ${workdir}/align

r1=${workdir}/orig/${filebase%bam}sorted.r1.fq
r2=${workdir}/orig/${filebase%bam}sorted.r2.fq
echo "$stardir/STAR --genomeDir $genome --readFilesIn $r1 $r2 --runThreadN 4 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 70000000000 --readFilesCommand cat --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --sjdbGTFfile $annotation --limitSjdbInsertNsj 2000000 --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --outSAMattrRGline ID::${tid}.${aid} --twopassMode Basic --outSAMmultNmax 1 > ${logdir}/${tid}.${aid}.align.log 2>&1 && touch $donefile_aln"
$stardir/STAR --genomeDir $genome --readFilesIn $r1 $r2 --runThreadN 4 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 70000000000 --readFilesCommand cat --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --sjdbGTFfile $annotation --limitSjdbInsertNsj 2000000 --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --outSAMattrRGline ID::${tid}.${aid} --twopassMode Basic --outSAMmultNmax 1 > ${logdir}/${tid}.${aid}.align.log 2>&1 && touch $donefile_aln

if [ ! -f $donefile_aln ]
then
    echo "Alignment step incomplete"
    rm -r ${workdir}
    rm ${result_dir}/running.lock
    exit 1
fi
mv ${workdir}/align/Aligned.sortedByCoord.out.bam $outfile_aln
mv ${workdir}/align/SJ.out.tab $outfile_jct

$samtools index $outfile_aln
echo "$(date) -- Alignment done"


###  STATS ###
##############
echo "$(date) -- Running alignment statistics in background"
python $stats_script $outfile_aln &


### EXPRESSION ###
##################
echo "$(date) -- Running expression counting in background"
outfile_cnt1="${outfile_aln%bam}ex_cnt.tab"
donefile_cnt1="${outfile_aln%bam}ex_cnt.done"
python $count_script -m -B -v -A ${outfile_aln} -a $annotation -o ${outfile_cnt1} > ${logdir}/${tid}.${aid}.exp_count.log 2>&1 && touch ${donefile_cnt1} &
outfile_cnt2="${outfile_aln%bam}ex_cnt.non_alt.tab"
donefile_cnt2="${outfile_aln%bam}ex_cnt.non_alt.done"
python $count_script -M -m -B -v -A ${outfile_aln} -a $annotation -o ${outfile_cnt2} > ${logdir}/${tid}.${aid}.exp_count.non_alt.log 2>&1 && touch ${donefile_cnt2} &

### SPLICING ###
################
echo "$(date) -- Running alternative splicing analysis"
mkdir -p ${workdir}/alt_splice
#mkdir -p ${workdir}/alt_splice/count
donefile_as="${workdir}/alt_splice/${tid}.${aid}.done"
echo "python ${as_script} -b ${outfile_aln} -o ${workdir}/alt_splice -a ${annotation_spladder} -v y -c 2 -M single -n 50 -P y -p n -T n --sparse_bam y --parallel $threads > ${logdir}/${tid}.${aid}.spladder.log && touch $donefile_as"
python ${as_script} -b ${outfile_aln} -o ${workdir}/alt_splice -a ${annotation_spladder} -v y -c 2 -M single -n 50 -P y -p n -T n --sparse_bam y --parallel $threads > ${logdir}/${tid}.${aid}.spladder.log 2>&1 && touch $donefile_as
if [ ! -f $donefile_as ]
then
    echo "SplAdder step incomplete"
    rm -r ${workdir}
    rm ${result_dir}/running.lock
    exit 1
fi
echo "$(date) -- Alternative splicing analysis done"

### EPILOGUE ###
################
echo "$(date) -- Waiting for all background processes to finish"
wait
echo "$(date) -- Moving relevant data to result directory"
mv ${workdir}/fastqc $result_dir/
mv ${workdir}/alt_splice $result_dir/
mv ${workdir}/align/*.hdf5 $result_dir/
mv ${workdir}/align/*.junctions $result_dir/
mv ${workdir}/align/*.stats $result_dir/
mv ${workdir}/align/*.tab $result_dir/
mv ${workdir}/align/*.done $result_dir/
touch $result_dir/pipeline.done
echo "$(date) -- Cleaning up"
rm -r ${workdir}
rm ${result_dir}/running.lock
echo "$(date) -- DONE"

