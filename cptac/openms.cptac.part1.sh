#!/bin/bash
usage="Usage: $0 -S sample_name -I fasta -C cptac_spectra_list -O output_dir"


base_dir=/cluster/work/grlab/projects/TCGA/PanCanAtlas/peptides_neoantigen/analysis_pancan

module load jdk/8u92
module load openms/2.1.0
msgfplus="/cluster/work/grlab/projects/TCGA/PanCanAtlas/peptides_neoantigen/sw/MSGFPlus/MSGFPlus.jar"
msgf_memory=80000
msgf_threads=6

cptac_ini="$base_dir/data/MSGF_iTRAQ.ini"




while getopts S:I:C:O:M: flag; do
    case $flag in
    S)
    sample=$OPTARG
    ;;
    I)
    fasta=$OPTARG
    ;;
    C)
    cptac_spectra_list=$OPTARG
    ;;
    O)
    output_dir=$OPTARG
    ;;
    M)
    cptac_ini=$OPTARG
    ;;
    \?)
    echo $usage >& 2
    exit -1
    ;;
    esac
done
shift $(( OPTIND - 1));



if [[ -z $sample ]]; then
    echo "Need to pass sample name to ${0}. Exiting ... "
    echo $usage
    exit -1
fi

if [[ -z $fasta ]]; then
    echo "Need to pass input fasta to ${0}. Exiting ... "
    echo $usage
    exit -1
fi

if [[ ! -s $fasta ]]; then
    echo "Input fasta $fasta does not exist. Exiting ... "
    echo $usage
    exit -1
fi

if [[ -z $cptac_spectra_list ]]; then
    echo "Need to pass file with all CPTAC spectra files (.mzML.gz) to ${0}. Exiting ... "
    echo $usage
    exit -1
fi

if [[ ! -s $cptac_spectra_list ]]; then
    echo "CPTAC spectra list file $cptac_spectra_list does not exist. Exiting ... "
    echo $usage
    exit -1
fi



if [[ -z $output_dir ]]; then
    echo "Need to pass output directory to ${0}. Exiting ..."
    echo $usage
    exit -1
fi

if [[ ! -d $output_dir ]]; then
    echo "Generating output directory $output_dir ..."
    mkdir -p $output_dir
fi

cptac_spectra_gz=`head -n$LSB_JOBINDEX $cptac_spectra_list | tail -n1`

if [[ ! -s $cptac_spectra_gz ]]; then
    echo "CPTAC spectra file $cptac_spectra_gz does not exist. Exiting ... "
    echo $usage
    exit -1
fi



me=`whoami`
master_time_start=`date +%s`
echo ">>>>>>>>>" >&2
echo Start OpenMS part I: `hostname --short` `date` >&2  

tmp_dir=/scratch/$me/raetsch_neoantigens_2016/cptac/${sample}/$RANDOM

mkdir -p $tmp_dir

bn=$(basename $cptac_spectra_gz)
label=$sample.${bn/.mzML.gz/}
output_prefix=$output_dir/$label
idXML=${output_prefix}.idXML
pi_idXML=${output_prefix}_pi.idXML


if [[ -s $pi_idXML ]] && [[ $pi_dXML -nt $idXML ]]; then
    echo File $pi_idXML exists and seems up-to-date. Skipping OpenMS partI ... 
    exit 1
fi

if [[ ! -s $idXML ]]; then
    echo "Unzip spectra file ... "
    cptac_spectra=$tmp_dir/${label}.mzML
    echo $cptac_spectra
    zcat $cptac_spectra_gz > $cptac_spectra
    ls -l $cptac_spectra

    echo "MSGFPlus" 
    echo MSGFPlusAdapter -ini $cptac_ini -in $cptac_spectra -out $idXML -database $fasta -executable $msgfplus -java_memory $msgf_memory -threads $msgf_threads
    MSGFPlusAdapter -ini $cptac_ini -in $cptac_spectra -out $idXML -database $fasta -executable $msgfplus -java_memory $msgf_memory -threads $msgf_threads
else
    echo File $idXML exists. Skipping MSGFPlusAdapter..
fi

if [[ ! -s $pi_idXML ]] || [[ $idXML -nt $pi_idXML ]]; then
    echo "PeptideIndexer"
    #PeptideIndexer -in $idXML -fasta $fasta -out $pi_idXML -allow_unmatched 
    PeptideIndexer -in $idXML -fasta $fasta -out $pi_idXML -allow_unmatched -enzyme:specificity 'semi'
else
    echo File $pi_idXML exists and seems up-to-date. Skipping PeptideIndexer ... 
fi

master_time_end=`date +%s`
(master_time_exec=`expr $(( $master_time_end - $master_time_start ))` 
echo "Finished running OpenMS MSGFPlusAdapter and PeptideIndexer on a chunk of $sample" >&2
echo "OpenMS CPTAC analysis part I completed in $master_time_exec seconds") >&2
echo End: `hostname --short` `date` >&2
echo "<<<<<<<<<" >&2

#rm -r $tmp_dir

