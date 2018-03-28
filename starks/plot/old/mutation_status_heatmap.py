import os
import sys
import pandas as pd
import numpy as np
import scipy as sp
import cPickle
import h5py

import sqtl
import utils
import config

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec

sys.path.append('/cluster/home/starks/git/tools-python/viz')
import heatmap

_LABELED_ENSG_PATH = '/cluster/work/grlab/projects/TCGA/PanCancer/rerun_plots/starks_sqtl/manhattan_sqtl_events_thold_80_v3_labeled_ensg.npy'
_FULL_META_PATH = '/cluster/work/grlab/projects/TCGA/PanCancer/metadata/2015-08-07/full_metadata.tsv'
_ENSG_TRANSLATE_PATH   = '/cluster/work/grlab/projects/ICGC/annotation/gencode.v19.annotation.hs37d5_chr.gtf.lookup.pickle'
_ENSG_TO_GENE = cPickle.load(open(_ENSG_TRANSLATE_PATH, 'r'))['ensembl2genename']
_VERSION = 3
_CISDATA_PATH = '/cluster/work/grlab/projects/TCGA/PanCancer/QTL_Analysis_results/analysis_10kset/2dmanhattan_information_v3/cisdataresults_test_v3.tsv'

def load_cis_data():
    cis_df = pd.read_csv(_CISDATA_PATH, sep='\t')
    if cis_df.columns[0].startswith('#'):
        cols = cis_df.columns.tolist()
        cols[0] = cols[0].replace('# ', '')
        cis_df.columns = cols
    cis_df = cis_df.set_index('gene-id(event)')
    return cis_df

def display_colname(colname, use_ensg=False):
    ensg, chrm, pos = colname.split('_')
    name_elems = list()
    ensg_split = ensg.split('-')
    gene_tags = list()
    for ensg_sng in ensg_split:
        if use_ensg:
            gene_tags.append(ensg_sng.split('.')[0])
        else:
            gene_tags.append(_ENSG_TO_GENE[ensg_sng.split('.')[0]])
    #ensg_name = '&'.join(gene_tags)
    ensg_name = gene_tags[0] # hot fix take first gene
    return '%s chr %s %s' %(ensg_name, chrm, pos)

def load_metadata():
    meta = pd.read_csv(_FULL_META_PATH, sep='\t')[['study', 'tcga_id']].set_index('tcga_id')
    meta = meta.drop([x for x in meta.index if not x.startswith('TCGA')])
    meta.index = meta.index.map(lambda x: _meta_index(x))
    return meta.iloc[:, 0].to_dict()

def _meta_index(_id):
    return '-'.join(_id.split('-')[:3])

def calc_mean_mutated(cnc_mut_status):
    mut_counts = np.sum(cnc_mut_status == 1)
    ref_counts = np.sum(cnc_mut_status == 2)
    return mut_counts.astype(float) / (mut_counts + ref_counts)

def _sort_colname_keys(col):
    '''Sort strings of '%s chr $chr $pos' by $chr then $pos'''
    chrm, pos = col.split()[-2:]
    return (int(chrm), int(pos))

def load_cisdata_mut_status():
    cis_df = load_cis_data()
    coords = cis_df.loc[:, ['snp_chrm', 'snp_pos', 'gene-id(mutation)']]
    coords.index = np.arange(coords.shape[0])
    coords = coords.drop_duplicates()
    coords = coords.dropna()
    coords.columns = ['chr', 'pos', 'gene']
    gt_file = h5py.File(sqtl._GENOTYPE_H5_PATH, 'r')
    gtdf = sqtl.get_mutation_info(gt_file, coords)
    return gtdf

def load_labeled_genes_mut_status():
    ensg_labeled = np.load(_LABELED_ENSG_PATH)
    gtdf = sqtl.load_gt_data()
    gtdf_cols = list()
    for ensg in ensg_labeled:
        mask = np.array(gtdf.columns.map(lambda x: ensg in x).tolist())
        assert mask.sum() > 0
        gtdf_cols.extend(gtdf.columns[mask].tolist())
    gtdf_cols = np.unique(gtdf_cols)
    mut_status = gtdf.loc[:, gtdf_cols]
    return mut_status

def plot_mut_status(mut_status, version, large_fig=False, ordering=None):
    meta = load_metadata()
    cnc_index = mut_status.index.map(lambda x: meta.get(_meta_index(x), 'UNK'))
    unq_cnc = np.unique(cnc_index)
    mut_mean_df = pd.DataFrame(index=unq_cnc, columns=mut_status.columns)

    for cnc in unq_cnc:
        cnc_mask = cnc_index == cnc
        mut_mean_df.loc[cnc] = calc_mean_mutated(mut_status.loc[cnc_mask])
    mut_mean_df = mut_mean_df.loc[:, mut_mean_df.mean(0) <= .5]

    if ordering is None:
        chrm = [int(x.split('_')[-2]) for x in mut_mean_df.columns]
        pos = [int(x.split('_')[-1]) for x in mut_mean_df.columns]
        order = np.lexsort((pos, chrm))
    elif ordering == 'cluster':
        _, order, _, _ = heatmap.cluster(mut_mean_df.values, dim1=False)
    elif ordering == 'entropy':
        order = np.argsort(sp.stats.entropy(mut_mean_df.values.astype(float)))
    elif ordering == 'maximum':
        order = np.argsort(mut_mean_df.values.astype(float).max(0))[::-1]
    else:
        raise ValueError('Unknown ordering argument')

    mut_mean_df = mut_mean_df.iloc[:, order]
    original_cols = mut_mean_df.columns
    ensg_cols = [display_colname(x, use_ensg=True) for x in mut_mean_df.columns]
    uniprot_cols = [display_colname(x, use_ensg=False) for x in mut_mean_df.columns]

    for col_type in ['ensg', 'uni']:
        fn_base = os.path.join('/cluster/work/grlab/projects/TCGA/PanCancer/rerun_plots', 'starks_sqtl', '%s_mutation_status_heatmaps'%version)
        if col_type == 'ensg':
            mut_mean_df.columns = ensg_cols
            fn_base = fn_base + '_ensg'
        else:
            mut_mean_df.columns = uniprot_cols
        fn_out = fn_base + '_ver%d.png'%_VERSION

        if large_fig:
            fig, ax = plt.subplots(figsize=(10, 50))
            cax = ax.imshow(mut_mean_df.values.astype(float).T, cmap='coolwarm', aspect='auto', vmax=1)
            ax.set_xticks(range(mut_mean_df.index.size))
            ax.set_xticklabels(mut_mean_df.index, rotation=90)
            ax.set_yticks(range(mut_mean_df.columns.size))
            ax.set_yticklabels(mut_mean_df.columns, rotation=0)

        else:
            fontsize = 14
            fig, ax = plt.subplots(figsize=(12, 8))
            cax = ax.imshow(mut_mean_df.values.astype(float).T, cmap='coolwarm', aspect='auto', vmax=1)
            ax.set_xticks(range(mut_mean_df.index.size))
            ax.set_xticklabels(mut_mean_df.index, rotation=90, fontsize=fontsize)
            ax.set_yticks(range(mut_mean_df.columns.size))
            ax.set_yticklabels(mut_mean_df.columns, fontsize=fontsize)
        cbar = fig.colorbar(cax, ticks=[0, .2, .4, .6, .8, 1])
        cbar.ax.set_title('% Mutated')
        plt.tight_layout()

        if not os.path.exists(os.path.dirname(fn_out)): os.makedirs(os.path.dirname(fn_out))
        print "Writing to %s" %fn_out
        plt.savefig(fn_out, bboxinches='tight', dpi=300)
        print "Writing to %s" %fn_out.replace('.png', '.pdf')
        plt.savefig(fn_out.replace('.png', '.pdf'), bboxinches='tight')
        plt.close()
    return

if __name__ == '__main__':
    do_labeled_genes = True
    do_cisdata = False

    if do_labeled_genes:
        labeled_genes_mut_status = load_labeled_genes_mut_status()
        plot_mut_status(labeled_genes_mut_status, 'labeled_genes')

    if do_cisdata:
        cisdata_mut_status = load_cisdata_mut_status()
        plot_mut_status(cisdata_mut_status, 'cisdata', large_fig=True, ordering='maximum')
        plot_mut_status(cisdata_mut_status, 'cisdata_entropy_ordering', large_fig=True, ordering='entropy')
