import os
import sys
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

def _print_embed_header(outdir, pca_ncomponents, tsne_lr_list, tsne_pp_list, do_tsne):
    if outdir is not None: print "Writing embeddings to %s" %outdir
    print "PCA embedding to %d dims"%pca_ncomponents
    if do_tsne:
        print "tSNE parameters:"
        print "  learning rates: %s" %str(tsne_lr_list)
        print "  perplexities: %s" %str(tsne_pp_list)

def run_embed_proc(data, index, outdir=None,
                   pca_ncomponents=50,
                   tsne_lr_list=[100, 500, 1000],
                   tsne_pp_list=[20, 40, 50, 75, 100, 200, 500],
                   do_tsne=True):

    assert data.shape[0] == index.size
    _print_embed_header(outdir, pca_ncomponents, tsne_lr_list, tsne_pp_list, do_tsne)
    pca_embeds_out = os.path.join(outdir, 'pca_embeds.tsv')
    pca_model_out = os.path.join(outdir, 'pca_model.npy')
    print "Starting PCA"
    pca = PCA(n_components=pca_ncomponents)
    pca_embeds = pca.fit_transform(data)
    assert pca_embeds.shape == (index.size, pca_ncomponents)
    print "Explained variance: %d%%" %(pca.explained_variance_ratio_.sum() * 100)

    if not outdir is None:
        pd.DataFrame(pca_embeds, index=index).to_csv(pca_embeds_out, sep='\t')
        np.save(pca_model_out, pca)
        print "Saving pca embedding to %s" %pca_embeds_out

    if do_tsne:
        if outdir is None: tsne_embed_dict = dict()
        for lr in tsne_lr_list:
            for pp in tsne_pp_list:
                print "\tFitting TSNE(perplexity=%d, learning_rate=%d)" %(pp, lr)
                tsne = TSNE(perplexity=pp, learning_rate=lr, n_iter=1000)
                tsne_embeds = tsne.fit_transform(pca_embeds)
                assert tsne_embeds.shape == (index.size, 2)

                if outdir is None:
                    tsne_embed_dict[(lr,pp)] = tsne_embeds
                else:
                    tsne_embeds_out = os.path.join(outdir, 'tsne_embeds_pp_%d_lr_%d.tsv'%(pp, lr))
                    tsne_model_out = os.path.join(outdir, 'tsne_model_pp_%d_lr_%d.npy'%(pp, lr))
                    tsne_embeds_df = pd.DataFrame(tsne_embeds, index=index).to_csv(tsne_embeds_out, sep='\t')
                    np.save(tsne_model_out, tsne)
                    print "Saving tsne embedding to %s" %tsne_embeds_out
    if outdir is None: return pca_embeds, tsne_embed_dict
    return
