#!/bin/bash
usage="Usage: $0 -S sample_name -I input_dir -O output_dir"

module load openms/2.1.0
ref="REFERENCE"

while getopts S:I:O: flag; do
    case $flag in
    S)
    sample=$OPTARG
    ;;
    I)
    input_dir=$OPTARG
    ;;
    O)
    output_dir=$OPTARG
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

if [[ -z $input_dir ]]; then
    echo "Need to pass input fasta to ${0}. Exiting ... "
    echo $usage
    exit -1
fi

if [[ ! -d $input_dir ]]; then
    echo "Input directory $input_dir does not exist. Exiting ... "
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



master_time_start=`date +%s`
echo ">>>>>>>>>" >&2
echo Start OpenMS: `hostname --short` `date` >&2  

echo "Merging input with background if given"
echo "Generating input fasta"
openms_input=$output_dir/$sample.openms_input.fa

cat $input_dir/${sample}*.fa $input_dir/${ref}*.fa > $openms_input

echo "Generating decoy database ..."
td_fasta=$output_dir/${sample}_td.fasta

if [[ ! -s $td_fasta ]]; then
    DecoyDatabase -in $openms_input -out $td_fasta
else
    echo "Decoy database $td_fasta exists. Skipping generation of decoy database."
fi

master_time_end=`date +%s`
(master_time_exec=`expr $(( $master_time_end - $master_time_start ))` 
echo "Finished generating decoy database on $sample" >&2
echo "OpenMS DecoyDatabases generation completed in $master_time_exec seconds") >&2
echo End: `hostname --short` `date` >&2
echo "<<<<<<<<<" >&2


