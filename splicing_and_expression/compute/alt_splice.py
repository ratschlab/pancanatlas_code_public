import os
import sys
import h5py
import numpy as np
import scipy as sp
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import config
import utils

DO_PCA = False
DO_TSNE = True and DO_PCA

DO_INDIVIDUAL_EVENTS = True
DO_COMBINED_EVENTS = False

PCA_NCOMPONENTS = 100
TSNE_LR_LIST = [1000]
TSNE_PP_LIST = [10, 50, 200, 500]

WHITELIST = utils.load_whitelist(config.whitelist_path)

FILTER_GTIDS = True
FILTER_GENE_PROCODE = True # only use protein coding genes

PLOT_DIR = os.path.join(config.plot_dir, 'altsplice')
EMBED_DIR = os.path.join(config.embed_dir, 'altsplice')

def run_preproc_tests(psi, psi_is_fin=None):
    '''Prints a table showing how many rows/cols left at each thold.
    '''
    if psi_is_fin is None: psi_is_fin = np.isfinite(psi)
    rows_pct_fin = np.mean(psi_is_fin, 1)
    assert rows_pct_fin.size == psi_is_fin.shape[0]
    cols_pct_fin = np.mean(psi_is_fin, 0)
    assert cols_pct_fin.size == psi_is_fin.shape[1]
    percentile_range = np.arange(0, 105, 5)
    for thold in np.linspace(0, 1, 11):
        n_cols_gt_thold = np.sum(cols_pct_fin >= thold)
        n_rows_gt_thold = np.sum(rows_pct_fin >= thold)
        pct_cols_gt_thold = np.mean(cols_pct_fin >= thold) * 100
        pct_rows_gt_thold = np.mean(rows_pct_fin >= thold) * 100
        filler = (thold, n_rows_gt_thold, pct_rows_gt_thold, n_cols_gt_thold, pct_cols_gt_thold)
        print "pct_finite >= %.2f  Rows left %6d (%2d%%) Columns left  %6d (%2d%%)" %filler
    print ""
    return

def replace_nans(psi, col_thold=.7, row_thold=.9):
    '''Replace nan entries in each row with mean of column

    Rows are samples, cols are events.
    '''
    psi_is_fin = np.isfinite(psi)
    print "Stats before filtering"
    run_preproc_tests(psi, psi_is_fin)

    print "Removing columns with < %d%% finite values" %int(100 * col_thold)
    cols_pct_fin = np.mean(psi_is_fin, 0)
    assert cols_pct_fin.size == psi_is_fin.shape[1]
    col_mask = cols_pct_fin > col_thold
    psi = psi[:, col_mask]
    psi_is_fin = psi_is_fin[:, col_mask]
    run_preproc_tests(psi, psi_is_fin)
    print ""

    print "Removing rows with < %d%% finite values" %int(100 * row_thold)
    rows_pct_fin = np.mean(psi_is_fin, 1)
    row_mask = rows_pct_fin > row_thold
    psi = psi[row_mask]
    psi_is_fin = psi_is_fin[row_mask]
    run_preproc_tests(psi, psi_is_fin)
    print ""

    print "Replacing nans with average of column"
    col_means = sp.nanmean(psi, 0)
    nan_ridx, nan_cidx = np.where(np.invert(psi_is_fin))
    psi[nan_ridx, nan_cidx] = col_means[nan_cidx]
    return psi, row_mask, col_mask

def preproc(psi, strains, gene_idx, conf_idx=None, ret_col_idx=False):
    '''Mask gene_idx to those in conf_idx.
    '''
    assert psi.shape == (strains.size, gene_idx.size)
    if not conf_idx is None:
        gene_idx = gene_idx[conf_idx]
        psi = psi[:, conf_idx]
    psi, row_mask, col_mask = replace_nans(psi)
    if ret_col_idx is not None:
        if conf_idx is None:
            col_idx = np.where(col_mask)[0]
        else:
            col_idx = conf_idx[col_mask]
    strains = strains[row_mask]
    gene_idx = gene_idx[col_mask]
    assert psi.shape == (strains.size, gene_idx.size)
    if ret_col_idx:
        return psi, strains, gene_idx, col_idx
    return psi, strains, gene_idx

def load_data(path, ret_col_idx=False):
    data = h5py.File(path, 'r')
    psi = data['psi'][:]
    print "\t psi.ngbytes = %.1f" %(psi.nbytes * 1e-9)
    gene_idx = data['gene_idx'][:]
    conf_idx = data['conf_idx'][:]
    strains = data['strains'][:]
    strains = clean_strain(strains)
    assert psi.shape[0] == strains.size
    data.close()
    print

    if WHITELIST is not None:
        mask = np.in1d(strains, WHITELIST)
        psi = psi[mask]
        strains = strains[mask]

    ret = preproc(psi, strains, gene_idx, conf_idx, ret_col_idx=ret_col_idx)
    return ret

def clean_strain(strain):
    strain = np.array([x.replace('.aligned', '').replace('.npz', '') for x in strain])
    return strain

def load_combined_events(map_event_to_file, thold=None):
    map_event_to_gene_index = dict()
    map_event_to_strains = dict()
    combined_index = None
    for key, path in map_event_to_file.items():
        _ , strains , gene_index = load_data(path)
        strains = clean_strain(strains)
        print "%s index size: %d" %(key, strains.size)
        combined_index = strains if combined_index is None else np.intersect1d(combined_index, strains)
        map_event_to_gene_index[key] = gene_index
        map_event_to_strains[key] = strains
    print "Combined index size: %d" %combined_index.size
    print ""

    combined_psi = np.zeros((combined_index.size, sum(v.size for v in map_event_to_gene_index.values())))
    combined_gidx = list()
    for event, path in map_event_to_file.items():
        print "Loading %s data from %s" %(event, path)
        psi, strains, gene_idx = load_data(path)
        sys.exit()
        strains = clean_strain(strains)
        assert np.all(np.in1d(combined_index, strains))
        assert psi.shape[0] == strains.size
        assert np.array_equal(gene_idx, map_event_to_gene_index[event])
        psi = pd.DataFrame(psi, index=strains).loc[combined_index].values
        combined_psi[:, len(combined_gidx):len(combined_gidx) + psi.shape[1]] = psi
        combined_gidx.extend(gene_idx)
        psi = None

    return combined_psi, combined_index, combined_gidx


if __name__ == '__main__':

    # load data
    map_event_to_file = {'exon_skip': config.alt_splice_exon_skip_path,
                         'intron_retention': config.alt_splice_intron_retention_path,
                         'alt_3prime': config.alt_splce_alt_3prime_path,
                         'alt_5prime': config.alt_splce_alt_5prime_path}

    total_gb = 0
    for key, value in map_event_to_file.items():
        shape = h5py.File(value, 'r')['psi'].shape
        gb = np.prod(shape) * 8 * 1e-9
        total_gb = total_gb + gb
        print "%s: %.1f Gb" %(key, gb)
    print "Total: %.1f Gb"%total_gb
    print ""

    if DO_INDIVIDUAL_EVENTS:
        for event, path in map_event_to_file.items():
            event_embed_dir = os.path.join(EMBED_DIR, event)
            if not os.path.exists(event_embed_dir): os.makedirs(event_embed_dir)

            input_cache = os.path.join(event_embed_dir, 'data.tsv')
            if os.path.exists(input_cache):
                print("Loading from cache: %s"%input_cache)
                df = pd.read_csv(input_cache, sep='\t', index_col=0)
                psi = df.values
                strains = df.index
                gene_idx = df.columns
            else:
                print "Loading %s data from %s" %(event, path)
                psi, strains, gene_idx = load_data(path)
                df = pd.DataFrame(psi, index=strains, columns=gene_idx)
                print("Caching to %s"%input_cache)
                df.to_csv(input_cache, sep='\t')

            if DO_PCA:
                print "PCA embedding to %d dims"%PCA_NCOMPONENTS
                pca = PCA(n_components=PCA_NCOMPONENTS)
                pca_embeds_norm = pca.fit_transform(psi)
                assert pca_embeds_norm.shape == (strains.size, PCA_NCOMPONENTS)
                pca_embeds_df = pd.DataFrame(pca_embeds_norm, index=strains)
                print "Explained variance: %d%%" %(pca.explained_variance_ratio_.sum() * 100)

                pca_embeds_out = os.path.join(event_embed_dir, 'pca_embeds.tsv')
                pca_model_out = os.path.join(event_embed_dir, 'pca_model.npy')
                if not os.path.exists(event_embed_dir): os.makedirs(event_embed_dir)
                pca_embeds_df.to_csv(pca_embeds_out, sep='\t')
                np.save(pca_model_out, pca)
                print "Saving pca embedding to %s" %pca_embeds_out

            if DO_TSNE:
                for lr in TSNE_LR_LIST:
                    for pp in TSNE_PP_LIST:
                        print "\tFitting TSNE(perplexity=%d, learning_rate=%d)" %(pp, lr)
                        tsne = TSNE(perplexity=pp, learning_rate=lr, n_iter=5000, random_state=4444)
                        tsne_embeds_norm = tsne.fit_transform(pca_embeds_df.values)
                        assert tsne_embeds_norm.shape == (strains.size, 2)
                        tsne_embeds_df = pd.DataFrame(tsne_embeds_norm, index=pca_embeds_df.index)

                        tsne_embeds_out = os.path.join(event_embed_dir, 'tsne_embeds_pp_%d_lr_%d.tsv'%(pp, lr))
                        tsne_model_out = os.path.join(event_embed_dir, 'tsne_model_pp_%d_lr_%d.npy'%(pp, lr))
                        tsne_embeds_df.to_csv(tsne_embeds_out, sep='\t')
                        np.save(tsne_model_out, tsne)
                        print "Saving tsne embedding to %s" %tsne_embeds_out
            psi = None; strains = None; gene_idx = None # for memory

    if DO_COMBINED_EVENTS:
        combined_psi, strains, _ = load_combined_events(map_event_to_file)
        event_embed_dir = os.path.join(EMBED_DIR, 'concatenated')
        if not os.path.exists(event_embed_dir): os.makedirs(event_embed_dir)

        input_cache = os.path.join(event_embed_dir, 'data.tsv')
        if os.path.exists(input_cache):
            print("Loading from cache: %s"%input_cache)
            df = pd.read_csv(input_cache, sep='\t', index_col=0)
            psi = df.values
            strains = df.index
            gene_idx = df.columns
        else:
            print "Loading %s data from %s" %(event, path)
            psi, strains, gene_idx = load_data(path)
            df = pd.DataFrame(psi, index=strains, columns=gene_idx)
            df.to_csv(input_cache, sep='\t')


        assert np.all(np.isfinite(combined_psi))
        if DO_PCA:
            print "PCA embedding to %d dims"%PCA_NCOMPONENTS
            pca = PCA(n_components=PCA_NCOMPONENTS)
            pca_embeds_norm = pca.fit_transform(combined_psi)
            combined_psi=None
            assert pca_embeds_norm.shape == (strains.size, PCA_NCOMPONENTS)
            pca_embeds_df = pd.DataFrame(pca_embeds_norm, index=strains)
            print "Explained variance: %d%%" %(pca.explained_variance_ratio_.sum() * 100)

            pca_embeds_out = os.path.join(event_embed_dir, 'pca_embeds.tsv')
            pca_model_out = os.path.join(event_embed_dir, 'pca_model.npy')
            pca_embeds_df.to_csv(pca_embeds_out, sep='\t')
            np.save(pca_model_out, pca)
            print "Saving pca embedding to %s" %pca_embeds_out

            if DO_TSNE:
                for lr in TSNE_LR_LIST:
                    for pp in TSNE_PP_LIST:
                        print "\tFitting TSNE(perplexity=%d, learning_rate=%d)" %(pp, lr)
                        tsne = TSNE(perplexity=pp, learning_rate=lr, n_iter=5000)
                        tsne_embeds_norm = tsne.fit_transform(pca_embeds_df.values)
                        assert tsne_embeds_norm.shape == (strains.size, 2)
                        tsne_embeds_df = pd.DataFrame(tsne_embeds_norm, index=pca_embeds_df.index)

                        tsne_embeds_out = os.path.join(event_embed_dir, 'tsne_embeds_pp_%d_lr_%d.tsv'%(pp, lr))
                        tsne_model_out = os.path.join(event_embed_dir, 'tsne_model_pp_%d_lr_%d.npy'%(pp, lr))
                        tsne_embeds_df.to_csv(tsne_embeds_out, sep='\t')
                        np.save(tsne_model_out, tsne)
                        print "Saving tsne embedding to %s" %tsne_embeds_out

