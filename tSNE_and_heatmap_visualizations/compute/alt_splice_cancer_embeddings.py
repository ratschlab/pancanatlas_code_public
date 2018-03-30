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

import alt_splice as preproc
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import config
import utils

DO_TSNE = True

DO_INDIVIDUAL_EVENTS = False
DO_COMBINED_EVENTS = True

PCA_NCOMPONENTS = 100
TSNE_LR_LIST = [100, 500, 1000]
TSNE_PP_LIST = [10, 20, 50, 100, 200, 500]

PLOT_DIR = os.path.join(config.plot_dir, 'altsplice')
EMBED_DIR = os.path.join(config.embed_dir, 'altsplice')


def run_pca_embeddings(psi, strains):
    print "PCA embedding to %d dims"%PCA_NCOMPONENTS
    pca = PCA(n_components=PCA_NCOMPONENTS)
    pca_embeds_norm = pca.fit_transform(psi)
    assert pca_embeds_norm.shape == (strains.size, min(psi.shape[0], PCA_NCOMPONENTS))
    pca_embeds_df = pd.DataFrame(pca_embeds_norm, index=strains)
    print "Explained variance: %d%%" %(pca.explained_variance_ratio_.sum() * 100)
    return pca_embeds_df, pca

def save_pca_embeddings(outdir, embed_df, pca_model):
    pca_embeds_out = os.path.join(outdir, 'pca_embeds.tsv')
    pca_model_out = os.path.join(outdir, 'pca_model.npy')
    if not os.path.exists(outdir): os.makedirs(outdir)
    pca_embeds_df.to_csv(pca_embeds_out, sep='\t')
    np.save(pca_model_out, pca_model)
    print "Saving pca embedding to %s" %pca_embeds_out
    return

def run_tsne_embeddings(pca_embeds_df, pp, lr):
    print "\tFitting TSNE(perplexity=%d, learning_rate=%d)" %(pp, lr)
    tsne = TSNE(perplexity=pp, learning_rate=lr, n_iter=5000)
    tsne_embeds_norm = tsne.fit_transform(pca_embeds_df.values)
    assert tsne_embeds_norm.shape == (pca_embeds_df.shape[0], 2)
    tsne_embeds_df = pd.DataFrame(tsne_embeds_norm, index=pca_embeds_df.index)
    return tsne_embeds_df, tsne

def save_tsne_embeddings(outdir, tsne_embeds_df, tsne, pp, lr):
    tsne_embeds_out = os.path.join(outdir, 'tsne_embeds_pp_%d_lr_%d.tsv'%(pp, lr))
    tsne_model_out = os.path.join(outdir, 'tsne_model_pp_%d_lr_%d.npy'%(pp, lr))
    tsne_embeds_df.to_csv(tsne_embeds_out, sep='\t')
    np.save(tsne_model_out, tsne)
    print "Saving tsne embedding to %s" %tsne_embeds_out
    return

def subset_to_cnc(psi, all_strains, cnc_strains):
    mask = np.in1d(all_strains, cnc_strains)
    if mask.sum() == 0: return None, None
    return psi[mask], all_strains[mask]

if __name__ == '__main__':

    # load data
    map_event_to_file = {'exon_skip': config.alt_splice_exon_skip_path,
                         'intron_retention': config.alt_splice_intron_retention_path,
                         'alt_3prime': config.alt_splce_alt_3prime_path,
                         'alt_5prime': config.alt_splce_alt_5prime_path}


    if DO_INDIVIDUAL_EVENTS:
        for event, path in map_event_to_file.items():
            print "Loading %s data from %s" %(event, path)
            psi, strains, gene_idx = preproc.load_data(path)
            metadata_df = utils.load_metadata_df(config.metadata_path, strains)
            cnc_groups = metadata_df.groupby('cnc')
            for cnc, cnc_md in cnc_groups:
                cnc_psi, cnc_strains = subset_to_cnc(psi, strains, cnc_md.index)
                if cnc_psi is None:
                    print "WARNING: Did not find %s samples in %s" %(cnc, event)
                cnc_event_embed_dir = os.path.join(EMBED_DIR, event, 'cancer_only', cnc)
                if not os.path.exists(cnc_event_embed_dir): os.makedirs(cnc_event_embed_dir)
                pca_embeds_df, pca = run_pca_embeddings(cnc_psi, cnc_strains)
                save_pca_embeddings(cnc_event_embed_dir, pca_embeds_df, pca)
                if DO_TSNE:
                    for lr in TSNE_LR_LIST:
                        for pp in TSNE_PP_LIST:
                            if pp > pca_embeds_df.shape[0]: continue
                            tsne_df, tsne = run_tsne_embeddings(pca_embeds_df, pp, lr)
                            assert tsne_df.index.equals(pca_embeds_df.index)
                            save_tsne_embeddings(cnc_event_embed_dir, tsne_df, tsne, pp, lr)
                cnc_psi = None; cnc_strains = None; gene_idx = None # for memory

    if DO_COMBINED_EVENTS:
        psi, strains, _ = preproc.load_combined_events(map_event_to_file)
        assert np.all(np.isfinite(psi))
        metadata_df = utils.load_metadata_df(config.metadata_path, strains)
        for cnc, cnc_md in metadata_df.groupby('cnc'):
            cnc_psi, cnc_strains = subset_to_cnc(psi, strains, cnc_md.index)
            if cnc_psi is None:
                print "WARNING: Did not find %s samples in combined" %cn
            cnc_event_embed_dir = os.path.join(EMBED_DIR, 'concatenated', 'cancer_only', cnc)
            if not os.path.exists(cnc_event_embed_dir): os.makedirs(cnc_event_embed_dir)
            pca_embeds_df, pca = run_pca_embeddings(cnc_psi, cnc_strains)
            save_pca_embeddings(cnc_event_embed_dir, pca_embeds_df, pca)
            if DO_TSNE:
                for lr in TSNE_LR_LIST:
                    for pp in TSNE_PP_LIST:
                        if pp > pca_embeds_df.shape[0]: continue
                        tsne_df, tsne = run_tsne_embeddings(pca_embeds_df, pp, lr)
                        assert tsne_df.index.equals(pca_embeds_df.index)
                        save_tsne_embeddings(cnc_event_embed_dir, tsne_df, tsne, pp, lr)


