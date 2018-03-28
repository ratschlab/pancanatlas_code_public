import os
import sys
import h5py
import argparse
import numpy as np
import pandas as pd

import matplotlib; matplotlib.use('Agg')
import seaborn as sns
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import utils

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'compute'))
import preproc


def name_outpath(max_events, only_pc, method, metric, collapsed):
    basename = '%d_high_var_events_heatmap.png' %max_events
    if collapsed: basename = 'collapsed_cnc_to_median_' + basename
    if only_pc: basename = 'protein_coding_' + basename
    return os.path.join(method + '_' + metric, basename)

def load_raw_counts(path):
    print "Loading expression counts (raw) from %s" %path
    expression_file = h5py.File(path, 'r')
    gids = expression_file['gids'][:]
    sids = expression_file['sids'][:]

    # counts are genes x gtids
    counts = expression_file['counts'][:]
    expression_file.close()
    return counts, gids, sids

def filter_counts(path, only_pc, max_events=None):
    counts, gids, sids = load_raw_counts(path)
    sids_mask = preproc.filter_sample_ids(sids)
    gids_mask = preproc.filter_genes(gids,
                                         counts=counts,
                                         only_protein_coding=only_pc,
                                         avg_lbound=FILTER_GENE_AVG_LBOUND)

    if not np.all(sids_mask):
        sids = sids[sids_mask]
        counts = counts[:, sids_mask]

    if not np.all(gids_mask):
        gids = gids[gids_mask]
        counts = counts[gids_mask]

    if max_events is not None:
        counts, gids = filter_to_high_var(counts.T, gids, max_events)

    # normalize
    counts = np.minimum(counts, np.percentile(counts, 99, axis=0))
    pct = np.percentile(counts, 75, axis=0)
    pct = np.maximum(1, pct)
    counts = counts / np.tile(pct, (counts.shape[0], 1))
    counts = np.log(counts + 1)
    return counts, gids, sids

def filter_to_high_var(data, columns, nkeep):
    '''Filter to the top nkeep high variance events
    '''
    if nkeep is None: return data, columns
    if nkeep <=1: nkeep = int(data.shape[1] * nkeep)
    var = np.var(data, axis=0)
    assert var.size == data.shape[1]
    keep_cols = np.argsort(var)[-nkeep:]
    data = data[:, keep_cols]
    columns = columns[keep_cols]
    return data, columns

def load_high_var_events_df(path, max_events=None, prefix=None, reset=None, only_pc=False):
    if prefix is None: cache_path = os.path.join(CACHE_DIR, 'counts_df.tsv')
    else: cache_path = os.path.join(CACHE_DIR, '%s_counts_df.tsv'%prefix)
    if reset is None: reset = RESET_DF_CACHE
    if reset or not os.path.exists(cache_path):
        print "Calculating counts"
        counts, gids, sids = filter_counts(path, only_pc=only_pc, max_events=max_events)
        df = pd.DataFrame(counts, index=sids, columns=gids)
        df.to_csv(cache_path, sep='\t')
    else:
        print "Reading counts from %s" %cache_path
        df = pd.read_csv(cache_path, sep='\t', index_col=0)
    return df

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
        print "Clustering df = " + str(df.shape)
        print "Calculating row linkage"
        row_linkage = hc.linkage(sp.distance.pdist(df.values), method=method, metric=metric)
        print "Calculating col linkage"
        col_linkage = hc.linkage(sp.distance.pdist(df.values.T), method=method, metric=metric)

        np.save(idx_linkage_cache_path, df.index)
        np.save(row_linkage_cache_path, row_linkage)
        np.save(col_linkage_cache_path, col_linkage)

    assert df.shape == (row_linkage.shape[0]+1, col_linkage.shape[0]+1)
    return row_linkage, col_linkage

def add_legend(graph, cnc_lut):
    '''Adds legend for row colors of graph
    '''
    for label, color in sorted(cnc_lut.items()):
        graph.ax_heatmap.bar(0, 0, color=color, label=label, linewidth=0)
    graph.ax_heatmap.legend(loc="upper center", ncol=4, bbox_to_anchor=(0.5, -.1))
    return

def collapse_cnc_to_median(combined_df, metadata_df):
    medians = dict()
    for cnc, cnc_df in metadata_df.groupby('cnc'):
        medians[cnc_df.index[0]] = combined_df.loc[cnc_df.index].median()
    return pd.DataFrame(medians).T

CACHE_DIR = os.path.expanduser(os.path.join('~/cache', 'rerun2018', 'expression_heatmap'))
if not os.path.exists(CACHE_DIR): os.makedirs(CACHE_DIR)

def define_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max_events', help='limit to number of events', type=int, default=5000)
    parser.add_argument('--only_protein_coding', help='subset to protein coding genes', action='store_true')
    parser.add_argument('--method', help='h-clustering method', default='ward')
    parser.add_argument('--metric', help='h-clustering metric', default='cosine')
    parser.add_argument('--reset_df_cache', help='clears cached counts', action='store_true')
    parser.add_argument('--reset_linkage_cache', help='clears cached clustering', action='store_true')
    parser.add_argument('--collapse_cnc_to_median', help='represent each cancer type as the median over samples', action='store_true')
    return parser

DEBUG = True

FILTER_GTIDS = True
FILTER_GENE_PROCODE = True # only use protein coding genes
FILTER_GENE_AVG_LBOUND = 50 # filter genes with avg count > value

if __name__ == '__main__':
    parser = define_parser()
    args = parser.parse_args()
    print args.__dict__
    # Global arguments
    MAX_EVENTS = args.max_events
    METHOD = args.method
    METRIC = args.metric
    RESET_DF_CACHE = args.reset_df_cache
    RESET_LINKAGE_CACHE = args.reset_linkage_cache
    ONLY_PC = args.only_protein_coding
    COLLAPSE_CNC_TO_MEDIAN = args.collapse_cnc_to_median

    if RESET_DF_CACHE: print "[WARNING]: Resetting df cache"
    if RESET_LINKAGE_CACHE: print "[WARNING]: Resetting linkage cache"

    PLOT_DIR = os.path.join(config.plot_dir, 'expression', 'heatmap')
    if not os.path.exists(PLOT_DIR): os.makedirs(PLOT_DIR)

    filtering_desc = '%d_high_var_events' %(MAX_EVENTS)
    if ONLY_PC: filtering_desc = 'protein_coding_%s' %filtering_desc
    title = filtering_desc + '_' + METHOD + '_' + METRIC
    if COLLAPSE_CNC_TO_MEDIAN: title = title + '_' + 'collapsed_cnc_to_median'

    # Load data, filter to MAX_EVENTS highest variance events (across all event_types)
    counts_df = load_high_var_events_df(config.expression_count_path, max_events=MAX_EVENTS, prefix=filtering_desc, only_pc=ONLY_PC, reset=RESET_DF_CACHE)
    metadata_df = utils.load_metadata_df(config.metadata_path, counts_df.index)

    # Get color scheme
    cnc_to_color = utils.load_color_scheme(config.color_scheme_path)
    cancer_types = cnc_to_color.keys()

    # Filter to samples within color scheme 
    index_subset = metadata_df.index[metadata_df['cnc'].isin(cancer_types)]
    counts_df = counts_df.loc[index_subset]
    metadata_df = metadata_df.loc[index_subset]
    if COLLAPSE_CNC_TO_MEDIAN:
        counts_df = collapse_cnc_to_median(counts_df, metadata_df)
    row_colors = metadata_df['cnc'].loc[counts_df.index].map(cnc_to_color)
    event_types = np.vectorize(lambda x: x.split('-')[0])(counts_df.columns) # event_type-idnum

    assert np.all(row_colors.index == counts_df.index)

    # Call clustering
    row_linkage, col_linkage = get_linkage(counts_df, method=METHOD, desc=title)

    # And finally plot the data
    sys.setrecursionlimit(100000)
    print "Plotting data ... "
    graph = sns.clustermap(counts_df,
                       row_colors=row_colors,
                       row_linkage=row_linkage, col_linkage=col_linkage,
                       cmap='BrBG')
    graph.ax_heatmap.axis('off')
    graph.ax_col_dendrogram.set_title("Expression %s Clustering" %title.replace('_', ' ').title())
    graph.ax_heatmap.set_xlabel("Events")
    graph.ax_heatmap.set_ylabel("Samples")
    graph.cax.set_title("psi")
    add_legend(graph, cnc_to_color)

    outpath_base = name_outpath(max_events=MAX_EVENTS, only_pc=ONLY_PC, method=METHOD, metric=METRIC, collapsed=COLLAPSE_CNC_TO_MEDIAN)
    outpath = os.path.join(PLOT_DIR, outpath_base)
    if not os.path.exists(os.path.dirname(outpath)): os.makedirs(os.path.dirname(outpath))
    print "Saving heatmap to: %s" %outpath
    plt.savefig(outpath, bbox_inches='tight')

