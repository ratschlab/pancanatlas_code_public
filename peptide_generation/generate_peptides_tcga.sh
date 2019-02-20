#!/usr/bin/env bash

export N_CORES=1

for i in $(seq 0 62)
do
    bsub -R "rusage[mem=64000]" -o ./cluster_logs/tcga_${i}_RESULT.txt -n $N_CORES -r -W 5:00 -J tcgaimmuno_${i} \
	python2 -u ./compute_peptide_lists.py --parallel_rank ${i} --batch_mode --run_id ccell_tcga
done
