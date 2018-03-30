import os
import sys
import h5py
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

import preproc
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import config
import utils

DO_CNC_BARS = False
DO_BOXPLOTS = False
DO_PCA = True
DO_TSNE = True and DO_PCA

PCA_NCOMPONENTS = 100
TSNE_LR_LIST = [500]
TSNE_PP_LIST = [10, 50, 200, 500]

FILTER_GTIDS = True
FILTER_GENE_PROCODE = True # only use protein coding genes
FILTER_GENE_AVG_LBOUND = 5 # filter genes with avg count > value

def load_data(processed_expression_count_path):
    # load expression counts
    print "Loading expression counts from %s" %processed_expression_count_path
    proc_expression_file = h5py.File(processed_expression_count_path, 'r')
    gtids = proc_expression_file['gtids'][:]
    gene_ids = proc_expression_file['expression']['expression']['gene_ids'][:]

    # counts are genes x gtids
    counts_norm = proc_expression_file['expression']['expression']['counts_norm'][:]
    proc_expression_file.close()
    return counts_norm, gene_ids, gtids

def load_raw_counts(expression_count_path):
    print "Loading expression counts (raw) from %s" %expression_count_path
    expression_file = h5py.File(expression_count_path, 'r')
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

def plot_cnc_bars(outpath, cnc):
    if not os.path.exists(PLOT_DIR): os.makedirs(PLOT_DIR)
    cnc_label, cnc_count = np.unique(cnc, return_counts=True)
    bar_pos = np.arange(cnc_label.size)
    fig, ax = plt.subplots()
    ax.bar(bar_pos, cnc_count, align='center', alpha=.5)
    ax.set_xticks(bar_pos)
    ax.set_xticklabels(cnc_label, rotation=90)
    ax.set_title('Samples per cancer type')
    print 'Writing %s' %outpath
    plt.savefig(outpath)
    return

def plot_box_raw(outpath, counts_raw):
    sample_indices = np.random.choice(counts_raw.shape[0], replace=False, size=25)
    fig, ax = plt.subplots()
    ax.boxplot(counts_raw[sample_indices].T)
    ax.set_ylabel('Log_10(count)')
    ax.set_xticks([])
    ax.set_title('Raw Counts')
    print "Writing %s" %outpath
    plt.savefig(outpath)
    return

def plot_box_norm(counts_norm):
    fig, ax = plt.subplots()
    counts_norm_sample = counts_norm[:, sample_indices]
    counts_norm_sample = counts_norm_sample + counts_norm_sample[counts_norm_sample > 0].min()
    ax.boxplot(-np.log10(counts_norm_sample))
    ax.set_ylabel('-log_10(norm_count)')
    ax.set_xticks([])
    ax.set_title('Normed Counts')
    boxplot_norm_path = os.path.join(PLOT_DIR, 'boxplot_norm_counts.png')
    print "Writing %s" %boxplot_norm_path
    plt.savefig(boxplot_norm_path)
    print "[scp call] scp hex:%s/* .; open *.png" %PLOT_DIR
    return

def compute_pca(counts_df):
    print "PCA embedding to %d dims"%PCA_NCOMPONENTS
    pca = PCA(n_components=PCA_NCOMPONENTS)
    pca_embeds = pca.fit_transform(counts_df.values)
    pca_embeds_df = pd.DataFrame(pca_embeds, index=counts_df.index)
    print "Explained variance: %d%%" %(pca.explained_variance_ratio_.sum() * 100)
    return pca_embeds_df, pca

def compute_tsne(embed_dir, pca_embeds_df):
    for lr in TSNE_LR_LIST:
        for pp in TSNE_PP_LIST:
            print "\tFitting TSNE(perplexity=%d, learning_rate=%d)" %(pp, lr)
            tsne = TSNE(perplexity=pp, learning_rate=lr, n_iter=5000, random_state=4444)
            tsne_embeds_norm = tsne.fit_transform(pca_embeds_df.values)
            tsne_embeds_df = pd.DataFrame(tsne_embeds_norm, index=pca_embeds_df.index)

            tsne_embeds_out = os.path.join(embed_dir, 'tsne_embeds_norm_pp_%d_lr_%d.tsv'%(pp, lr))
            tsne_model_out = os.path.join(embed_dir, 'tsne_model_norm_pp_%d_lr_%d.npy'%(pp, lr))
            print "Saving tsne embedding to %s" %tsne_embeds_out
            tsne_embeds_df.to_csv(tsne_embeds_out, sep='\t')
            np.save(tsne_model_out, tsne)

def run_qc_plots(plot_dir, counts_df, metadata_df=None):
    if DO_CNC_BARS:
        assert metadata_df is not None
        plot_cnc_bars(cnc)

    if DO_BOXPLOTS:
        plot_box(counts_df)
    return

def run_embed_routines(counts_df):
    return

if __name__ == '__main__':
    PLOT_DIR = os.path.join(config.plot_dir, 'expression')
    EMBED_DIR = os.path.join(config.embed_dir, 'expression')
    if not os.path.exists(PLOT_DIR): os.makedirs(PLOT_DIR)
    if not os.path.exists(EMBED_DIR): os.makedirs(EMBED_DIR)

    input_cache = os.path.join(EMBED_DIR, 'data.tsv')
    if os.path.exists(input_cache):
        print("Loading from cache: %s"%input_cache)
        counts_raw_df = pd.read_csv(input_cache, sep='\t', index_col=0)
    else:
        print "Loading expression count data from %s" %config.expression_count_path
        wl = utils.load_whitelist(config.whitelist_path)
        counts_raw, gids, sids = load_raw_counts(config.expression_count_path)
        wl_mask = np.in1d(sids, wl)
        assert counts_raw.shape[1] == sids.size
        counts_raw = counts_raw[:, wl_mask]
        sids = sids[wl_mask]
        counts_raw, gids, sids = filter_counts(counts_raw, gids, sids)
        assert np.all(np.isfinite(counts_raw))
        counts_raw = np.log10(counts_raw + 1).T
        counts_raw_df = pd.DataFrame(counts_raw, index=sids, columns=gids)
        # df is sids x gids
        counts_raw_df.to_csv(input_cache, sep='\t')
    metadata_df = utils.load_metadata_df(config.metadata_path)
    metadata_df = utils.translate_tcga_to_strain_index(metadata_df, counts_raw_df.index)

    if DO_CNC_BARS:
        bar_outpath = os.path.join(PLOT_DIR, 'cancer_type_barplot.png')
        cnc = metadata_df['cnc'].loc[counts_raw_df.index].values
        plot_cnc_bars(bar_outpath, cnc)

    if DO_BOXPLOTS:
        box_outpath = os.path.join(PLOT_DIR, 'expression_boxplot.png')
        plot_box_raw(box_outpath, counts_raw_df.values)

    if DO_PCA:
        pca_embeds_raw_out = os.path.join(EMBED_DIR, 'pca_embeds.tsv')
        pca_model_raw_out = os.path.join(EMBED_DIR, 'pca_model.npy')

        pca_embeds_raw_df, pca_model_raw = compute_pca(counts_raw_df)
        print "Saving pca embedding to %s" %pca_embeds_raw_out
        pca_embeds_raw_df.to_csv(pca_embeds_raw_out, sep='\t')
        np.save(pca_model_raw_out, pca_model_raw)

    if DO_TSNE:
        assert DO_PCA
        compute_tsne(EMBED_DIR, pca_embeds_raw_df)

