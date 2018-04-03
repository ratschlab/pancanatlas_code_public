#!/bin/bash

#!/bin/bash
usage="Usage: $0 -S sample_name -O working_dir"

#msgfplus="/cluster/home/tonora/software/msgfplus/MSGFPlus.jar"
#msgf_memory=80000
#msgf_threads=6



while getopts S:O: flag; do
    case $flag in
    S)
    sample=$OPTARG
    ;;
    O)
    working_dir=$OPTARG
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

if [[ -z $working_dir ]]; then
    echo "Need to pass output directory to ${0}. Exiting ..."
    echo $usage
    exit -1
fi

if [[ ! -d $working_dir ]]; then
    echo "$working_dir does not exist. Exiting ..."
    exit -1
fi


input_dir=$working_dir/chunks

if [[ ! -d $input_dir ]]; then
    echo Expecting subdirectory chunks in working directory $working_dir. Cannot find it. Exiting ...
    exit -1
fi

module load openms/2.1.0
me=`whoami`
master_time_start=`date +%s`
echo ">>>>>>>>>" >&2
echo Start OpenMS part II: `hostname --short` `date` >&2  

echo "Merge pi.idXML" 
#f1=`ls ${input_dir}/*_[fF]1_pi.idXML`
#f10=`ls ${input_dir}/*[Ff]10_pi.idXML` #f10 or F10
f10=`ls ${input_dir}/*10_pi.idXML` #f10 or F10
#tmp_merged_file=${f10/_[Ff]10_pi/_pi}

helper=${f10%_*}
tmp_merged_file=${helper%_*}_pi.idXML
merged_file=${tmp_merged_file/chunks\//}
echo $tmp_merged_file
echo $merged_file


if [[ ! -s $merged_file ]] || [[ $f10 -nt $merged_file ]]; then
    echo "IDMerger -in $input_dir/*_[Ff]*[0-9]_pi.idXML -out $merged_file"
    IDMerger -in $input_dir/*_[Ff]*[0-9]_pi.idXML -out $merged_file
else
    echo File $merged_file exists. Skipping IDMerger..
fi

echo "FalseDiscoveryRate"
fd_idXML=${merged_file/_pi/_fdr}
if [[ ! -s $fd_idXML ]] || [[ $merged_file -nt $fd_idXML ]]; then
    echo "FalseDiscoveryRate -in $merged_file -out $fd_idXML"
    FalseDiscoveryRate -in $merged_file -out $fd_idXML
else
    echo File $fd_idXML exists and seems up-to-date. Skipping FalseDiscoveryRate ... 
fi

echo "IDFilter"
cutoff=0.05
h=${cutoff/\./}
filtered_idXML=${merged_file/_pi/_fdr$h}
if [[ ! -s $filtered_idXML ]] || [[ $fd_idXML -nt $filtered_idXML ]]; then
    echo "IDFilter -in $fd_idXML -out $filtered_idXML -score:pep $cutoff"
    IDFilter -in $fd_idXML -out $filtered_idXML -score:pep $cutoff
else
    echo File $filtered_idXML exists and seems up-to-date. Skipping IDFilter ...
fi

echo "FileInfo"
fileinfo_out=${filtered_idXML/idXML/info}
if [[ ! -s $fileinfo_out ]] || [[ $filtered_idXML -nt $fileinfo_out ]]; then
    echo "FileInfo -in $filtered_idXML | tee $fileinfo_out"
    FileInfo -in $filtered_idXML | tee $fileinfo_out
else
    echo File $fileinfo_out exists and seems up-to-date. Skipping FileInfo ... 
fi


echo "TextExporter"
filtered_csv=${filtered_idXML/idXML/csv}
if [[ ! -s $filtered_csv ]] || [[ $filtered_idXML -nt $filtered_csv ]]; then
    echo "TextExporter -in $filtered_idXML -out $filtered_csv"
    TextExporter -in $filtered_idXML -out $filtered_csv
else
    echo File $filtered_csv exists and seems up-to-date. Skipping TextExporter ...
fi


master_time_end=`date +%s`
(master_time_exec=`expr $(( $master_time_end - $master_time_start ))` 
echo "Finished running OpenMS IDMerger, FalseDiscoveryRate & IDFilter on $sample" >&2
echo "OpenMS CPTAC analysis part II completed in $master_time_exec seconds") >&2
echo End: `hostname --short` `date` >&2
echo "<<<<<<<<<" >&2


