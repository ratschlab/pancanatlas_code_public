import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import config
import utils

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plot_utils

import numpy as np
import pandas as pd

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc

RESET_CACHE = True

def get_linkage(df, desc, method='ward', metric='euclidean'):
    '''Heirarchical linkage on df rows/cols, or load from cache

        checks indicies after loading from cache
    '''
    sample_linkage_cache_path = os.path.join(config.CACHE_DIR, '%s_sample_linkage.npy'%desc)
    event_linkage_cache_path = os.path.join(config.CACHE_DIR, '%s_event_linkage.npy'%desc)
    idx_linkage_cache_path = os.path.join(config.CACHE_DIR, '%s_idx_linkage.npy'%desc)
    if RESET_CACHE and os.path.exists(sample_linkage_cache_path) and os.path.exists(event_linkage_cache_path):
        print "WARNING: Overwriting linkage from: \n  %s\n  %s" %(sample_linkage_cache_path, event_linkage_cache_path)
    if not RESET_CACHE and os.path.exists(sample_linkage_cache_path) and os.path.exists(event_linkage_cache_path):
        print "Reading linkage from: \n  %s\n  %s" %(sample_linkage_cache_path, event_linkage_cache_path)
        idx_linkage = np.load(idx_linkage_cache_path)
        assert np.all(df.index == idx_linkage)
        sample_linkage = np.load(sample_linkage_cache_path)
        event_linkage = np.load(event_linkage_cache_path)
    else:
        print "Clustering df = " + str(df.shape)
        print "Calculating row linkage"
        sample_linkage = hc.linkage(sp.distance.pdist(df.values, metric=metric), method=method, metric=metric)
        print "Calculating col linkage"
        event_linkage = hc.linkage(sp.distance.pdist(df.values.T, metric=metric), method=method, metric=metric)

        print "Writing linkage to: \n  %s\n  %s" %(sample_linkage_cache_path, event_linkage_cache_path)
        np.save(idx_linkage_cache_path, df.index)
        np.save(sample_linkage_cache_path, sample_linkage)
        np.save(event_linkage_cache_path, event_linkage)

    assert df.shape == (sample_linkage.shape[0]+1, event_linkage.shape[0]+1)
    return sample_linkage, event_linkage

def filter_to_high_var(data, columns, nkeep):
    '''Filter to the top nkeep high variance columns
    '''
    if nkeep is None: return data, columns
    if nkeep <=1: nkeep = int(data.shape[1] * nkeep)
    var = np.var(data, axis=0)
    assert var.size == data.shape[1]
    keep_cols = np.argsort(var)[-nkeep:]
    return keep_cols

def add_legend(graph, cnc_lut):
    '''Adds legend for row colors of graph
    '''
    for label, color in sorted(cnc_lut.items()):
        graph.ax_heatmap.bar(0, 0, color=color, label=label, linewidth=0)
    graph.ax_heatmap.legend(loc="upper center", ncol=4, bbox_to_anchor=(0.5, -.1))
    return

def collapse_to_median(df, meta):
    medians = dict()
    for label, subser in meta.groupby(meta):
        medians[subser.index[0]] = df.loc[subser.index].median()
    return pd.DataFrame(medians).T

def main(df, outdir, desc, color_loader, run_representative):
    '''Runs all tasks on a single embedding.

    embed_dir: location of pca & tsne embeddings
    plot_dir: directory to write plots
    desc: identifies embedding (used for e.g. plot titles)
    '''
    print("clustermap: %s" %desc)
    if not os.path.exists(outdir): os.makedirs(outdir)

    assert desc.lower().startswith('altsplice') or desc.lower().startswith('expression')
    is_altsplice = desc.lower().startswith('altsplice')
    if not is_altsplice:
        assert df.values.max() > 1, "99.999999% Sure this is not psi data"
        print("Clipping 99th percentile")
        df.iloc[:] = np.minimum(df.values, np.percentile(df.values, 99, axis=0))

    method = 'ward'
    metric = 'cosine'

    # get colors
    cat_series, color_lut = color_loader(df)
    colors = cat_series.map(color_lut)
    df = df.loc[colors.index]

    # filter to high var
    keep_cols = filter_to_high_var(df.values, df.columns, MAX_EVENTS)
    df = df.iloc[:, keep_cols]

    cluster_desc = '_'.join(['%d_high_var_events'%MAX_EVENTS, method, metric])
    sample_desc = desc.strip().replace(' ', '_').lower() + '_' + cluster_desc

    sample_linkage, event_linkage = get_linkage(df, sample_desc, method=method, metric=metric)
    outpath = os.path.join(outdir, cluster_desc + '_clustermap.png')
    plot_heatmap(outpath, df, sample_linkage, colors, event_linkage, desc, color_lut)

    if run_representative:
        medians = collapse_to_median(df, cat_series)
        rep_colors = colors.loc[medians.index]
        rep_cluster_desc = cluster_desc + '_reps'
        rep_desc = desc.strip().replace(' ', '_').lower() + '_' + rep_cluster_desc
        print("clustermap: %s" %desc)
        rep_sample_linkage, rep_event_linkage = get_linkage(medians, rep_desc, method=method, metric=metric)
        rep_outpath = os.path.join(outdir, rep_cluster_desc + '_clustermap.png')
        plot_heatmap(rep_outpath, medians, rep_sample_linkage, rep_colors, rep_event_linkage, rep_desc, color_lut)
    return

def plot_pairwise_dists(linkage, colors, color_lut, desc):
    import ipdb; ipdb.set_trace()
    return

def plot_heatmap(outpath, df, sample_linkage, sample_colors, event_linkage, desc, sample_color_lut):

    assert desc.lower().startswith('altsplice') or desc.lower().startswith('expression')
    is_altsplice = desc.lower().startswith('altsplice')

    sys.setrecursionlimit(100000)
    print "Plotting data ... "
    graph = sns.clustermap(df.T,
                       col_colors=sample_colors,
                       col_linkage=sample_linkage, row_linkage=event_linkage,
                       cmap = sns.cubehelix_palette(as_cmap=True))
    graph.ax_heatmap.axis('off')
    graph.ax_col_dendrogram.set_title("%s Clustering" %' '.join(desc.split('_')).title())
    graph.ax_heatmap.set_xlabel("Events")
    graph.ax_heatmap.set_ylabel("Samples")
    if is_altsplice: graph.cax.set_title("psi")
    else: graph.cax.set_title("log(counts)")
    add_legend(graph, sample_color_lut)
    plot_utils.save(outpath)
    return

def load_mrna_colors(df):
    key = 'Subtype_mRNA'
    metadata_df = utils.load_metadata_df(config.metadata_path, df.index)
    metadata_df, subtype_names = utils.append_subtype(config.subtype_path, metadata_df)
    assert key in subtype_names and key in metadata_df.columns

    subtype_series = metadata_df[key].loc[df.index].dropna().apply(lambda x: x.lower())
    value_counts = subtype_series.value_counts()
    subtype_classes = value_counts[value_counts > 30].index
    subtype_series = subtype_series[subtype_series.isin(subtype_classes)]
    colors = utils.load_color_scheme(config.color_scheme_path).values()
    color_lut = dict(zip(sorted(subtype_classes), colors))
    return subtype_series, color_lut

def load_cnc_colors(df):
    metadata_df = utils.load_metadata_df(config.metadata_path, df.index)

    # Get color scheme
    cnc_to_color = utils.load_color_scheme(config.color_scheme_path)
    cancer_types = cnc_to_color.keys()
    cnc_colors = metadata_df['cnc'].map(cnc_to_color)
    return metadata_df['cnc'], cnc_to_color

# Tasks
CNC = True
mRNA = False

# Filtering
MAX_EVENTS = 5000
DEBUG = False

if __name__ == '__main__':
    path_list = list()
    outdir_list = list()
    desc_list = list()

    # Add Expression
    if True:
        path_list.append(os.path.join(config.embed_dir, 'expression', 'data.tsv'))
        outdir_list.append(os.path.join(config.plot_dir, 'expression', 'heatmaps'))
        desc_list.append('Expression')

    # Add AltSplice
    if True:
#        altsplice_event_list= ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime']
        altsplice_event_list= ['exon_skip']
        for event in altsplice_event_list:
            path_list.append(os.path.join(config.embed_dir, 'altsplice', event, 'data.tsv'))
            outdir_list.append(os.path.join(config.plot_dir, 'altsplice', event, 'heatmap'))
            desc_list.append('AltSplice %s'%event.title())

    for i, (path, outdir, desc) in enumerate(zip(path_list, outdir_list, desc_list)):
        assert os.path.exists(path), path
        print("[%d] %s: %s"%(i, desc, path))
        print("[%d] %s Plots:  %s"%(i, desc, outdir))

    for i, (path, outdir, desc) in enumerate(zip(path_list, outdir_list, desc_list)):
        print "Loading %s" %path
        try:
            data = utils.load_large_df(path.replace('.tsv', ''))
        except IOError:
            data = pd.read_csv(path, sep='\t', index_col=0)

    if CNC:
        main(data, os.path.join(outdir, 'cnc'), desc + ' CNC', load_cnc_colors, run_representative=True)
    if mRNA:
        main(data, os.path.join(outdir, 'mrna'), desc + ' mRNA', load_mrna_colors, run_representative=True)

