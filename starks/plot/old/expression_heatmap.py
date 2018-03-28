import os
import sys
import numpy as np
import h5py
import pandas as pd

import matplotlib; matplotlib.use('Agg')
import seaborn as sns
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import config
import utils

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'compute'))
import preproc

FILTER_GTIDS = True
FILTER_GENE_PROCODE = True # only use protein coding genes
FILTER_GENE_AVG_LBOUND = 50 # filter genes with avg count > value

def load_raw_counts(config):
    print "Loading expression counts (raw) from %s" %config.expression_count_path
    expression_file = h5py.File(config.expression_count_path, 'r')
    gids = expression_file['gids'][:]
    sids = expression_file['sids'][:]

    # counts are genes x gtids
    counts_raw = expression_file['counts'][:]
    expression_file.close()
    return counts_raw, gids, sids

def filter_counts(counts, gene_ids, gtids):
    gtids_mask = preproc.filter_sample_ids(gtids)
    gene_ids_mask = preproc.filter_genes(gene_ids,
                                         counts=counts,
                                         only_protein_coding=FILTER_GENE_PROCODE,
                                         avg_lbound=FILTER_GENE_AVG_LBOUND)

    if not np.all(gtids_mask):
        gtids = gtids[gtids_mask]
        counts = counts[:, gtids_mask]

    if not np.all(gene_ids_mask):
        gene_ids = gene_ids[gene_ids_mask]
        counts = counts[gene_ids_mask]

    return counts, gene_ids, gtids

def load_data(config):
    counts_raw, gids, sids = load_raw_counts(config)
    counts_raw, gids, sids = filter_counts(counts_raw, gids, sids)
    assert np.all(np.isfinite(counts_raw))
    counts_raw = np.log10(counts_raw + 1).T

    if TEST:
        print "[TEST MODE]: subsampling counts uniformly"
        row_idx = np.random.choice(counts_raw.shape[0], size=counts_raw.shape[0]/10)
        col_idx = np.random.choice(counts_raw.shape[1], size=counts_raw.shape[1]/10)
        counts_raw = counts_raw[row_idx]
        counts_raw = counts_raw[:, col_idx]
        sids = sids[row_idx]
        gids = gids[col_idx]
    counts_raw_df = pd.DataFrame(counts_raw, index=sids, columns=gids)
    counts_raw_df = counts_raw_df.reset_index().drop_duplicates(subset='index', keep='last').set_index('index')
    return counts_raw_df

def load_cnc_colors(config, counts_raw_df, max_cnc=20):
    '''Assigns a color to each sample by its cancer type.
    '''
    metadata_df = utils.load_metadata_df(config.metadata_path)
    metadata_df = utils.translate_tcga_to_strain_index(metadata_df, counts_raw_df.index)
    metadata_df = metadata_df.loc[metadata_df.index & counts_raw_df.index]

    top_cnc = metadata_df['cnc'].value_counts()[:max_cnc].index
    metadata_df = metadata_df.loc[metadata_df['cnc'].isin(top_cnc)]
    counts_raw_df = counts_raw_df.loc[metadata_df.index]

    cnc_pal = sns.color_palette('tab20', max_cnc)
    cnc_pal = cnc_pal[::2] + cnc_pal[1::2]
    cnc_lut = dict(zip(top_cnc, cnc_pal))
    cnc_colors = metadata_df.loc[counts_raw_df.index]['cnc'].map(cnc_lut)
    return cnc_colors, counts_raw_df, zip(top_cnc, cnc_pal)

def get_linkage(counts_raw_df, **kwargs):
    row_linkage_cache_path = os.path.join(CACHE_DIR, '%s_row_linkage.npy'%METHOD)
    col_linkage_cache_path = os.path.join(CACHE_DIR, '%s_col_linkage.npy'%METHOD)

    if not TEST and CACHE and os.path.exists(row_linkage_cache_path) and os.path.exists(col_linkage_cache_path):
        "Reading from: \n  %s\n  %s" %(row_linkage_cache_path, col_linkage_cache_path)
        row_linkage = np.load(row_linkage_cache_path)
        col_linkage = np.load(col_linkage_cache_path)
    else:
        print "Calculating row linkage"
        row_linkage = hc.linkage(sp.distance.pdist(counts_raw_df.values), **kwargs)
        print "Calculating col linkage"
        col_linkage = hc.linkage(sp.distance.pdist(counts_raw_df.values.T), **kwargs)
        if not TEST:
            np.save(row_linkage_cache_path, row_linkage)
            np.save(col_linkage_cache_path, col_linkage)
    return row_linkage, col_linkage

def _name_outpath():
    plot_dir = os.path.join(config.plot_dir, 'raw_expression')
    if not os.path.exists(plot_dir): os.makedirs(plot_dir)
    if TEST:
        outpath = os.path.join(plot_dir, 'test_%s_heatmap.png'%METHOD)
    else:
        outpath = os.path.join(plot_dir, '%s_heatmap.png'%METHOD)
    return outpath

def plot(counts, **kwargs):
    g = sns.clustermap(counts, **kwargs)
    g.ax_col_dendrogram.set_title("Raw Expression %s Clustering" %METHOD.title())
    g.ax_heatmap.set_xlabel("Genes")
    g.ax_heatmap.set_ylabel("Samples")
    g.cax.set_title("log(counts)")
    return g

def add_legend(graph, cnc_lut):
    for label, color in cnc_lut:
        graph.ax_heatmap.bar(0, 0, color=color, label=label, linewidth=0)
    graph.ax_heatmap.legend(loc="upper center", ncol=4, bbox_to_anchor=(0.5, -.1))
    return g

METHOD = "ward"
CACHE = True
CACHE_DIR = os.path.expanduser(os.path.join('~/cache', 'rerun2018', 'expression_heatmap'))
if not os.path.exists(CACHE_DIR): os.makedirs(CACHE_DIR)
TEST = False

if __name__ == '__main__':
    print "Clustering method: %s" %METHOD

    counts_raw_df = load_data(config)
    cnc_colors, counts_raw_df, cnc_lut = load_cnc_colors(config, counts_raw_df)
    row_linkage, col_linkage = get_linkage(counts_raw_df, method=METHOD)
    assert counts_raw_df.shape == (row_linkage.shape[0]+1, col_linkage.shape[0]+1)
    outpath = _name_outpath()

    print "Plotting data ... "
    g = plot(counts_raw_df,
             row_colors=cnc_colors.values,
             row_linkage=row_linkage, col_linkage=col_linkage,
             xticklabels=False, yticklabels=False,
             linewidths=0)

    g = add_legend(g, cnc_lut)
    print "Saving heatmap to: %s" %outpath
    sns.plt.savefig(outpath, bbox_inches='tight')

