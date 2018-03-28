import os
import sys
import argparse

BASEDIR = os.path.dirname(os.path.dirname(__file__))
sys.path.append(BASEDIR)
import compute.alt_splice as preproc
import config
import utils

import numpy as np
import pandas as pd

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc

def load_cnc_colors(config, df, max_cnc=20):
    '''Assigns a color to each sample by its cancer type.
    '''
    metadata_df = utils.load_metadata_df(config.metadata_path)
    metadata_df = utils.translate_tcga_to_strain_index(metadata_df, df.index)
    joint_index = metadata_df.index & df.index
    if joint_index.size != df.index.size:
        print "[WARNING]: %d samples do not have metadata"%(joint_index.size - df.index.size)
    metadata_df = metadata_df.loc[joint_index]

    top_cnc = metadata_df['cnc'].value_counts()[:max_cnc].index
    metadata_df = metadata_df.loc[metadata_df['cnc'].isin(top_cnc)]
    df = df.loc[metadata_df.index]
    print "Dropping %d rare cancer type samples" %(joint_index.size - df.index.size)

    cnc_pal = sns.color_palette('tab20', max_cnc)
    cnc_pal = cnc_pal[::2] + cnc_pal[1::2]
    cnc_lut = dict(zip(top_cnc, cnc_pal))
    cnc_colors = metadata_df.loc[df.index]['cnc'].map(cnc_lut)
    return cnc_colors, df, zip(top_cnc, cnc_pal)

def map_col_colors(event_names, lut):
    col_colors = [lut[ev.split('-')[0]] for ev in event_names]
    col_colors = pd.Series(col_colors, index=event_names)
    return col_colors

def get_linkage(df, desc, method='ward', metric='euclidean', reset=None):
    '''Heirarchical linkage on df rows/cols, or load from cache

        checks indicies after loading from cache
    '''
    prefix = desc.strip().replace(' ', '_').lower()
    row_linkage_cache_path = os.path.join(CACHE_DIR, '%s_row_linkage.npy'%prefix)
    col_linkage_cache_path = os.path.join(CACHE_DIR, '%s_col_linkage.npy'%prefix)
    idx_linkage_cache_path = os.path.join(CACHE_DIR, '%s_idx_linkage.npy'%prefix)
    if reset is None: reset = RESET_LINKAGE_CACHE
    if not reset and os.path.exists(row_linkage_cache_path) and os.path.exists(col_linkage_cache_path):
        print "Reading linkage from: \n  %s\n  %s" %(row_linkage_cache_path, col_linkage_cache_path)
        idx_linkage = np.load(idx_linkage_cache_path)
        assert np.all(df.index == idx_linkage)
        row_linkage = np.load(row_linkage_cache_path)
        col_linkage = np.load(col_linkage_cache_path)
    else:
        print "Calculating row linkage"
        row_linkage = hc.linkage(sp.distance.pdist(df.values), method=method, metric=metric)
        print "Calculating col linkage"
        col_linkage = hc.linkage(sp.distance.pdist(df.values.T), method=method, metric=metric)

        np.save(idx_linkage_cache_path, df.index)
        np.save(row_linkage_cache_path, row_linkage)
        np.save(col_linkage_cache_path, col_linkage)

    assert df.shape == (row_linkage.shape[0]+1, col_linkage.shape[0]+1)
    return row_linkage, col_linkage

def add_legend(graph, lut):
    '''Adds legend for row colors of graph

    lut: dict mapping labels to colors
    graph: sns.clustermap object
    '''
    for label, color in sorted(lut.items()):
        graph.ax_heatmap.bar(0, 0, color=color, label=label, linewidth=0)
    graph.ax_heatmap.legend(loc="upper center", ncol=4, bbox_to_anchor=(0.5, -.1))
    return

def add_col_legend(graph, col_colors):
    for label, color in sorted(col_colors.items()):
        graph.ax_col_dendrogram.bar(0, 0, color=color, label=label, linewidth=0)
    graph.ax_col_dendrogram.legend(loc='center', ncol=len(col_colors))
    return

CACHE_DIR = os.path.expanduser('~/cache/alt_splice_heatmap')
if not os.path.exists(CACHE_DIR): os.makedirs(CACHE_DIR)

def define_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max_events', help='limit to number of events', type=int, default=5000)
    parser.add_argument('--only_protein_coding', help='subset to protein coding genes', action='store_true')
    parser.add_argument('--method', help='h-clustering method', default='ward')
    parser.add_argument('--metric', help='h-clustering metric', default='euclidean')
    parser.add_argument('--reset_df_cache', help='clears cached counts', action='store_true')
    parser.add_argument('--reset_linkage_cache', help='clears cached clustering', action='store_true')
    event_choices = ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime', 'combined']
    parser.add_argument('--event_type', help='pick which event to consider', choices=event_choices, default='combined')
    parser.add_argument('--collapse_cnc_to_median', help='represent each cancer type as the median over samples', action='store_true')
    return parser

DEBUG = True

if __name__ == '__main__':
    parser = define_parser()
    args = parser.parse_args()
    print args.__dict__
    # Global arguments
    METHOD = args.method
    METRIC = args.metric
    RESET_DF_CACHE = args.reset_df_cache
    RESET_LINKAGE_CACHE = args.reset_linkage_cache
    ONLY_PC = args.only_protein_coding
    COLLAPSE_CNC_TO_MEDIAN = args.collapse_cnc_to_median

    if RESET_DF_CACHE: print "[WARNING]: Resetting df cache"
    if RESET_LINKAGE_CACHE: print "[WARNING]: Resetting linkage cache"

    # Get color scheme
    row_lut = utils.load_color_scheme(config.color_scheme_path)
    cancer_types = cnc_to_color.keys()

    # Filter to samples within color scheme 
    index_subset = metadata_df.index[metadata_df['cnc'].isin(cancer_types)]
    combined_df = combined_df.loc[index_subset]
    metadata_df = metadata_df.loc[index_subset]


    # cluster rows & cols
    # make plot
    # return plotting object, linkages

    # input args:
    #   combined_df -- events x samples
    #   row_colors, col_colors
    #   row_lut, col_lut
    #   {method, metric} or linkages -- clustering args


    def make_heatmap(data_df,
                     row_colors=None, col_colors=None,
                     row_lut=None, col_lut=None,
                     method=None, metric=None,
                     row_linkage=None, col_linkage=None):
    '''Cluster map from data_df with row, col colors
    '''
    specd_linkages = not (row_linkags is None and col_linkage is None)
    specd_clustering = not (method is None and metric is None)
    assert specd_linkages or specd_clustering, "need to specify method & metric or one of row_linkage, col_linkage"

    if not specd_linkages:
        if specd_clustering: print("Warning: ignoring cluster specifications")
        row_linkage, col_linkage = get_linkage(data_df, method=method, metric=metric)

    # And finally plot the data
    sys.setrecursionlimit(100000)
    print "Plotting data ... "
    graph = sns.clustermap(data_df,
                       row_colors=row_colors, col_colors=col_colors,
                       row_linkage=row_linkage, col_linkage=col_linkage)
    graph.ax_heatmap.axis('off')
    if not row_lut is None: add_legend(graph, row_lut)
    return graph, row_linkage, col_linkage



altsplice_event_list= ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime']
# path, dtype, dtype_sub
if dtype == 'altsplice':
    data_df = aspp.load_high_var_events_single_df(path, max_events=5000, only_pc=only_pc)
else:
    data_df = expp.load_high_var_events_single_df(path, max_events=5000, only_pc=only_pc)

metadata_df = utils.load_metadata_df(config.metadata_path, data_df.index)
