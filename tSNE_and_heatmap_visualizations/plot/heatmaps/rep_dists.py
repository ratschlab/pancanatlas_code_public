import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import config
import utils

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plot_utils

import numpy as np
import pandas as pd
import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as sch
import scipy.stats as spst

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

def collapse_to_median(df, meta):
    medians = pd.DataFrame(index=meta.unique(), columns=df.columns, dtype=float)
    for label, subser in meta.groupby(meta):
        medians.loc[label] = df.loc[subser.index].median()
    return medians

def filter_to_high_var(data, columns, nkeep):
    '''Filter to the top nkeep high variance columns
    '''
    if nkeep is None: return data, columns
    if nkeep <=1: nkeep = int(data.shape[1] * nkeep)
    var = np.var(data, axis=0)
    assert var.size == data.shape[1]
    keep_cols = np.argsort(var)[-nkeep:]
    return keep_cols

def heatmap_dists_with_dendro(data, norm=False, labels=None, metric='euclidean', method='ward'):

    fig = plt.figure(figsize=(7 * 1.30, 7 * 1.25))
    gs = gridspec.GridSpec(ncols=3, nrows=2, height_ratios=[.25, 1], width_ratios=[.25, 1, .05], hspace=0)
    dend_top_ax = fig.add_subplot(gs[0,1])
    hmap_ax     = fig.add_subplot(gs[1,1])
    cbar_ax     = fig.add_subplot(gs[1,2])
    dend_top_ax.set_axis_off()

    if labels is None:
        try:
            labels = data.index
        except AttributeError:
            pass

    n = data.shape[0]
    assert labels is None or len(labels) == n

    dists = ssd.pdist(data, metric=metric)
    linkage = sch.linkage(dists, metric=metric, method=method)
    dendro = sch.dendrogram(linkage, ax=dend_top_ax, color_threshold=0, above_threshold_color='black')
    order = dendro['leaves']
    sq_form_dists = ssd.squareform(dists)[order][:, order]
    assert sq_form_dists.shape == (n,n)

    if norm:
        sq_form_dists = spst.zscore(sq_form_dists, axis=None)
        sq_form_dists *= -1
        cmap = plt.get_cmap('cubehelix')
        vmin = -4
        vmax = 4
    else:
        cmap = plt.get_cmap()
        vmin = None
        vmax = None
    hmap = hmap_ax.imshow(sq_form_dists, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
    hmap_ax.set_xticks(np.arange(n))
    hmap_ax.set_yticks(np.arange(n))
    if labels is not None:
        hmap_ax.set_xticklabels(labels[order], rotation=90)
        hmap_ax.set_yticklabels(labels[order])
    cb = plt.colorbar(hmap, cax=cbar_ax)
    return

def heatmap_dists(data, norm=False, labels=None, metric='euclidean', method='ward'):
    fig, (ax, cax) = plt.subplots(ncols=2,figsize=(7 * 1.05 ,7),
                                  gridspec_kw={"width_ratios":[1, 0.05]})

    if labels is None:
        try:
            labels = data.index
        except AttributeError:
            pass

    n = data.shape[0]
    assert labels is None or len(labels) == n

    dists = ssd.pdist(data, metric=metric)
    linkage = sch.linkage(dists, metric=metric, method=method)
    dendro = sch.dendrogram(linkage, no_plot=True)
    order = dendro['leaves']
    sq_form_dists = ssd.squareform(dists)[order][:, order]
    assert sq_form_dists.shape == (n,n)

    hmap = ax.imshow(sq_form_dists, aspect='auto')
    ax.set_xticks(np.arange(n))
    ax.set_yticks(np.arange(n))
    if labels is not None:
        ax.set_xticklabels(labels[order], rotation=90)
        ax.set_yticklabels(labels[order])
    cb = plt.colorbar(hmap, cax=cax)
    return fig, (ax, cax)


# Tasks
CNC = True
mRNA = False

# Filtering
MAX_EVENTS = 5000
DEBUG = False
NORM = True

if __name__ == '__main__':
    path_list = list()
    outdir_list = list()
    desc_list = list()

    # Add Expression
    if True:
        path = os.path.join(config.embed_dir, 'expression', 'data.tsv')
        outdir = os.path.join(config.plot_dir, 'expression', 'heatmaps')
        if not os.path.exists(outdir): os.makedirs(outdir)
        desc = 'Expression'
        if NORM: desc = 'Normalized ' + desc
        try:
            df = utils.load_large_df(path.replace('.tsv', ''))
        except IOError:
            df = pd.read_csv(path, sep='\t', index_col=0)

        df.iloc[:] = np.minimum(df.values, np.percentile(df.values, 99, axis=0))
        keep_cols = filter_to_high_var(df.values, df.columns, MAX_EVENTS)
        df = df.iloc[:, keep_cols]
        metadata_df = utils.load_metadata_df(config.metadata_path, df.index)
        medians = collapse_to_median(df, metadata_df['cnc'])
        heatmap_dists_with_dendro(medians, norm=NORM)
        outpath = os.path.join(outdir, desc.lower().replace(' ', '_') +'_rep_dists_heatmap.png')
        plot_utils.save(outpath, do_pdf=True)

    # Add AltSplice
    if False:
        altsplice_event_list= ['alt_3prime', 'alt_5prime', 'intron_retention', 'exon_skip']
        for event in altsplice_event_list:
            path = os.path.join(config.embed_dir, 'altsplice', event, 'data.tsv')
            outdir = os.path.join(config.plot_dir, 'altsplice', event, 'heatmap')
            if not os.path.exists(outdir): os.makedirs(outdir)
            desc = 'AltSplice %s'%event.title()
            if NORM: desc = 'Normalized ' + desc
            print desc
            print "Loading %s" %path
            try:
                df = utils.load_large_df(path.replace('.tsv', ''))
            except IOError:
                df = pd.read_csv(path, sep='\t', index_col=0)

            keep_cols = filter_to_high_var(df.values, df.columns, MAX_EVENTS)
            df = df.iloc[:, keep_cols]
            metadata_df = utils.load_metadata_df(config.metadata_path, df.index)
            medians = collapse_to_median(df, metadata_df['cnc'])
            heatmap_dists_with_dendro(medians, metric='cosine', norm=NORM)
            outpath = os.path.join(outdir, desc.lower().replace(' ', '_') +'_rep_dists_heatmap_TEST.png')
            plot_utils.save(outpath, do_pdf=True)

