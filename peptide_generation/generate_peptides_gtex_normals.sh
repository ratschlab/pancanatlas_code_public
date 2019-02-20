#!/usr/bin/env bash

export N_CORES=1

for i in $(seq 0 62)
do
    bsub -R "rusage[mem=64000]" -o ./cluster_logs/gtex_${i}_RESULT.txt -n $N_CORES -r -W 5:00 -J tcgaimmuno_gtex_${i} \
	python2 -u ./compute_peptide_lists.py --parallel_rank ${i} --batch_mode --gtex_normals --run_id ccell_gtex
done
