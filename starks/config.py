import os
import sys
import pandas as pd

#FIGURE DUMP: https://docs.google.com/presentation/d/1gP5RdMU2vUn4lhv0nyWJ9FXkpxjkC967c7CNNYFMywo/edit?usp=sharing

# TASKS: https://docs.google.com/document/d/1vdCF0s746wxyo2e2ivbP7RbVcCeocyf9dFQ7uj9pAEU/edit
PROJECT_DIR = '/cluster/work/grlab/projects/TCGA/PanCancer'
CACHE_DIR = os.path.join(PROJECT_DIR, 'cache')

color_scheme_path = os.path.join(PROJECT_DIR, 'annotation/color_schemes/PanCanAtlasTumors_color_coded_by_organ_system_20170302.tsv')

# TCGA barcode: https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode
metadata_path = os.path.join(PROJECT_DIR, 'metadata/2015-08-07/full_metadata.overlap_RNA_WXS.RNA.tsv')
exome_metadata_path = os.path.join(PROJECT_DIR, 'metadata/2015-08-07/full_metadata.overlap_RNA_WXS.WXS.tsv')
subtype_path = os.path.join(PROJECT_DIR, 'annotation/subtypes/TCGASubtype.20170308.clean.tsv')
staging_info = '/cluster/work/grlab/projects/TCGA/PanCancer/annotation/clinical/clinical_PANCAN_patient_with_followup.stage_only.tsv'
gtex_metadata_path = '/cluster/work/grlab/projects/ICGC/gtex_tables/SraRunTable.txt'
smoking_age_path = '/cluster/work/grlab/projects/TCGA/PanCancer/annotation/clinical/clinical_PANCAN_patient_with_followup.age_and_smoke.tsv'

# generate PCA / tSNE plots for
#expression_count_path = os.path.join(PROJECT_DIR, 'rerun_hdf5/expression_counts.hdf5')
#processed_expression_count_path = os.path.join(PROJECT_DIR, 'hdf5_10Kset/phenotypes.hdf5r10_s2000_V100.hdf5')
expression_count_path = os.path.join(PROJECT_DIR, 'rerun2018_hdf5', 'expression_counts.hdf5')

#alt_splice_dir = os.path.join(PROJECT_DIR, 'rerun_alt_splice')
alt_splice_dir = os.path.join(PROJECT_DIR, 'rerun2018_alt_splice')
alt_splice_exon_skip_path = os.path.join(alt_splice_dir, 'merge_graphs_exon_skip_C2.counts.hdf5')
alt_splice_intron_retention_path = os.path.join(alt_splice_dir, 'merge_graphs_intron_retention_C2.counts.hdf5')
alt_splce_alt_3prime_path = os.path.join(alt_splice_dir, 'merge_graphs_alt_3prime_C2.counts.hdf5')
alt_splce_alt_5prime_path = os.path.join(alt_splice_dir, 'merge_graphs_alt_5prime_C2.counts.hdf5')
alt_splice_alt_3prime_path = alt_splce_alt_3prime_path
alt_splice_alt_5prime_path = alt_splce_alt_5prime_path
psi_cache_dir = os.path.join(alt_splice_dir, 'psi_caches')

alt_splice_gtex_dir = os.path.join(PROJECT_DIR, 'rerun_alt_splice.GTEx')
alt_splice_exon_skip_gtex_path = os.path.join(alt_splice_gtex_dir, 'merge_graphs_exon_skip_C3.counts.hdf5')
alt_splice_intron_retention_gtex_path = os.path.join(alt_splice_gtex_dir, 'merge_graphs_intron_retention_C3.counts.hdf5')
alt_splce_alt_3prime_gtex_path = os.path.join(alt_splice_gtex_dir, 'merge_graphs_alt_3prime_C3.counts.hdf5')
alt_splce_alt_5prime_gtex_path = os.path.join(alt_splice_gtex_dir, 'merge_graphs_alt_5prime_C3.counts.hdf5')

genotype_somatic_path = os.path.join(PROJECT_DIR, 'hdf5_10Kset/somaticVariants.hdf5')
genotype_germline_path = os.path.join(PROJECT_DIR, 'hdf5_10Kset/mergedFiles_clean.hdf5')

plot_dir = os.path.join(PROJECT_DIR, 'rerun2018_plots')
embed_dir = os.path.join(PROJECT_DIR, 'rerun2018_embeds')

fn_protein_coding_list = os.path.join(PROJECT_DIR, 'annotation/gencode.v19.annotation.hs37d5_chr_prtnCoding_list.txt')
census_gene_path = os.path.join(PROJECT_DIR, 'annotation/gencodeV14.v7.pancan_subset.ensembleID.list')

libsize_path = os.path.join(PROJECT_DIR, 'rerun2018_hdf5/expression_counts.whitelisted.libsize.tsv')
rnadeg_path  = os.path.join(PROJECT_DIR, 'rerun2018_qc/degradation/deg_files.tsv')
whitelist_path = os.path.join(PROJECT_DIR, 'rerun2018_qc/sample_whitelist.txt')


def load_protein_coding_genes():
    data = pd.read_csv(fn_protein_coding_list, delimiter='\t', usecols=[0]).values.ravel()
    return data


