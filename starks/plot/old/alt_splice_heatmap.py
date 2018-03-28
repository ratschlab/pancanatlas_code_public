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

def load_data(path):
   psi, strains, gene_idx = preproc.load_data(path)
   gene_idx = gene_idx.astype(int)
   assert psi.shape == (strains.size, gene_idx.size)
   return pd.DataFrame(psi, index=strains, columns=gene_idx)

def name_outpath(max_events, only_pc, method, metric, collapsed):
    basename = '%d_high_var_events_heatmap.png' %max_events
    if collapsed: basename = 'collapsed_cnc_to_median_' + basename
    if only_pc: basename = 'protein_coding_' + basename
    return os.path.join(method + '_' + metric, basename)

def get_combined_df(map_event_to_file, max_events=None, prefix=None, index=None, reset=None, only_pc=False):
    '''Read and filter psi data or load from cache
    '''
    if prefix is None:
        cache_path = os.path.join(CACHE_DIR, 'combined_df.tsv')
    else:
        cache_path = os.path.join(CACHE_DIR, '%s_combined_df.tsv'%prefix)

    if reset is None: reset = RESET_DF_CACHE
    if reset or not os.path.exists(cache_path):
        print "Calculating psi"
        combined_df = load_high_var_events_combined_df(map_event_to_file, max_events, index=index)
        combined_df.to_csv(cache_path, sep='\t')
    else:
        print "Reading psi from %s" %cache_path
        combined_df = pd.read_csv(cache_path, sep='\t', index_col=0)
    return combined_df

def load_high_var_events_single_df(path, max_events=None, index=None, only_pc=False):
    '''Load data from each altsplice category, filter to high-var events.

        map_event_to_file: dict mapping event_type -> data csv path
        max_events: filter to top max_events highest variance events
        index: filter to events in the index
        only_pc: if true, subset genes to only protein coding
    '''

    psi, strains, gene_idx = preproc.load_data(path)
    strains = preproc.clean_strain(strains)
    if index is not None:
        mask = np.in1d(strains, index)
        assert mask.sum() > 0
        strains = strains[mask]
        psi = psi[mask]
    if only_pc:
        psi, gene_idx = _filter_to_protein_coding(psi, gene_idx)
    if max_events is not None:
        psi, gene_idx = filter_to_high_var(psi, gene_idx, max_events)
    df = pd.DataFrame(psi, index=strains, columns=gene_idx)
    return df

def load_high_var_events_combined_df(map_event_to_file, max_events=None, index=None, condense=True, only_pc=False):
    '''Load data from each altsplice category, filter to high-var events.

        map_event_to_file: dict mapping event_type -> data csv path
        max_events: filter to top max_events highest variance events
        index: filter to events in the index
        condense: if true, return top events across all event types
        only_pc: if true, subset genes to only protein coding
    '''
    df_list = list()
    for event, path in map_event_to_file.items():
        print "Loading %s psi from %s" %(event, path)
        df = load_high_var_events_single_df(path, max_events=max_events, index=index, only_pc=only_pc)
        df.columns = [event + '-' + str(int(x)) for x in df.columns]
        df_list.append(df)

    df = pd.concat(df_list, axis=1, join='inner')
    assert df.shape[0] > 0
    if max_events is not None:
        values = df.values
        columns = df.columns
        index = df.index
        df = None
        if condense: values, columns = filter_to_high_var(values, columns, max_events)
        df = pd.DataFrame(values, index=index, columns=columns)
    return df

def _filter_to_protein_coding(psi, gene_idx, gene_names):
    pc_ensgs = np.loadtxt(config.fn_protein_coding_list, delimiter='\t', dtype=str, usecols=[0])
    import ipdb; ipdb.set_trace()
    return psi, gene_idx


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

def add_legend(graph, cnc_lut):
    '''Adds legend for row colors of graph
    '''
    for label, color in sorted(cnc_lut.items()):
        graph.ax_heatmap.bar(0, 0, color=color, label=label, linewidth=0)
    graph.ax_heatmap.legend(loc="upper center", ncol=4, bbox_to_anchor=(0.5, -.1))
    return

def add_col_legend(graph, col_colors):
    for label, color in sorted(col_colors.items()):
        graph.ax_col_dendrogram.bar(0, 0, color=color, label=label, linewidth=0)
    graph.ax_col_dendrogram.legend(loc='center', ncol=len(col_colors))
    return

def process(combined_df, outpath=None):
    if combined_df.shape[1] > MAX_EVENTS:
        index = combined_df.index
        data, columns = filter_to_high_var(combined_df.values, combined_df.columns, MAX_EVENTS)
        combined_df = None; combined_df = pd.DataFrame(data, index=index, columns=columns)
    cnc_colors, combined_df, cnc_lut = load_cnc_colors(config, combined_df, MAX_CNC)
    print "psi shape: %s" %str(combined_df.shape)
    row_linkage, col_linkage = get_linkage(combined_df, method=METHOD)
    assert combined_df.shape == (row_linkage.shape[0]+1, col_linkage.shape[0]+1)

    sys.setrecursionlimit(100000)
    print "Plotting data ... "
    g = plot(combined_df,
             row_colors=cnc_colors.values,
             row_linkage=row_linkage, col_linkage=col_linkage,
             xticklabels=False, yticklabels=False,
             linewidths=0)

    g = add_legend(g, cnc_lut)
    return

def collapse_cnc_to_median(combined_df, metadata_df):
    medians = dict()
    for cnc, cnc_df in metadata_df.groupby('cnc'):
        medians[cnc_df.index[0]] = combined_df.loc[cnc_df.index].median()
    return pd.DataFrame(medians).T

def old_args():
    MAX_EVENTS = 1000
    METHOD = 'ward' # h-clustering method
    METRIC = 'cosine'
    RESET_DF_CACHE = False
    RESET_LINKAGE_CACHE = False

CACHE_DIR = os.path.expanduser('~/cache/alt_splice_heatmap')
if not os.path.exists(CACHE_DIR): os.makedirs(CACHE_DIR)


def define_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max_events', help='limit to number of events', type=int, default=5000)
    parser.add_argument('--only_protein_coding', help='subset to protein coding genes', action='store_true')
    parser.add_argument('--method', help='h-clustering method', default='ward')
    parser.add_argument('--metric', help='h-clustering metric', default='cosine')
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
    MAX_EVENTS = args.max_events
    METHOD = args.method
    METRIC = args.metric
    RESET_DF_CACHE = args.reset_df_cache
    RESET_LINKAGE_CACHE = args.reset_linkage_cache
    ONLY_PC = args.only_protein_coding
    COLLAPSE_CNC_TO_MEDIAN = args.collapse_cnc_to_median

    # Local arguments
    etype_desc = args.event_type

    if RESET_DF_CACHE: print "[WARNING]: Resetting df cache"
    if RESET_LINKAGE_CACHE: print "[WARNING]: Resetting linkage cache"

    # event_type to path of raw psi values
    map_etype_to_file = {'exon_skip': config.alt_splice_exon_skip_path,
                         'intron_retention': config.alt_splice_intron_retention_path,
                         'alt_3prime': config.alt_splce_alt_3prime_path,
                         'alt_5prime': config.alt_splce_alt_5prime_path}
    if etype_desc != 'concatenated': map_etype_to_file = {etype_desc: map_etype_to_file[etype_desc]}

    PLOT_DIR = os.path.join(config.plot_dir, 'altsplice', etype_desc, 'heatmap')
    if not os.path.exists(PLOT_DIR): os.makedirs(PLOT_DIR)

#    map_etype_to_file = {'intron_retention': config.alt_splice_intron_retention_path}
#    etype_desc = 'intron_retention'

    filtering_desc = '%s_%d_high_var_events' %(etype_desc, MAX_EVENTS)
    if ONLY_PC: filtering_desc = 'protein_coding_%s' %filtering_desc
    title = filtering_desc + '_' + METHOD + '_' + METRIC
    if COLLAPSE_CNC_TO_MEDIAN: title = title + '_' + 'collapsed_cnc_to_median'

    # Load data, filter to MAX_EVENTS highest variance events (across all event_types)
    combined_df = get_combined_df(map_etype_to_file, MAX_EVENTS, prefix=filtering_desc, only_pc=ONLY_PC, reset=RESET_DF_CACHE)
    metadata_df = utils.load_metadata_df(config.metadata_path, combined_df.index)

    # Get color scheme
    cnc_to_color = utils.load_color_scheme(config.color_scheme_path)
    cancer_types = cnc_to_color.keys()
    col_cmap = np.array(sns.color_palette('Set2', len(map_etype_to_file.keys())))
    col_cmap_lut = dict(zip(map_etype_to_file.keys(), col_cmap))

    # Filter to samples withing color scheme 
    index_subset = metadata_df.index[metadata_df['cnc'].isin(cancer_types)]
    combined_df = combined_df.loc[index_subset]
    metadata_df = metadata_df.loc[index_subset]
    if COLLAPSE_CNC_TO_MEDIAN:
        combined_df = collapse_cnc_to_median(combined_df, metadata_df)
    row_colors = metadata_df['cnc'].loc[combined_df.index].map(cnc_to_color)
    event_types = np.vectorize(lambda x: x.split('-')[0])(combined_df.columns) # event_type-idnum

    col_colors = map_col_colors(combined_df.columns, col_cmap_lut)
    assert np.all(row_colors.index == combined_df.index)

    # Call clustering
    row_linkage, col_linkage = get_linkage(combined_df, method=METHOD, desc=title)

    # And finally plot the data
    sys.setrecursionlimit(100000)
    print "Plotting data ... "
    graph = sns.clustermap(combined_df,
                       row_colors=row_colors, col_colors=col_colors,
                       row_linkage=row_linkage, col_linkage=col_linkage)
    graph.ax_heatmap.axis('off')
    graph.ax_col_dendrogram.set_title("AltSplice %s Clustering" %title.replace('_', ' ').title())
    graph.ax_heatmap.set_xlabel("Events")
    graph.ax_heatmap.set_ylabel("Samples")
    graph.cax.set_title("psi")
    add_legend(graph, cnc_to_color)
    add_col_legend(graph, col_cmap_lut)

    outpath_base = name_outpath(max_events=MAX_EVENTS, only_pc=ONLY_PC, method=METHOD, metric=METRIC, collapsed=COLLAPSE_CNC_TO_MEDIAN)
    outpath = os.path.join(PLOT_DIR, outpath_base)
    if not os.path.exists(os.path.dirname(outpath)): os.makedirs(os.path.dirname(outpath))
    print "Saving heatmap to: %s" %outpath
    plt.savefig(outpath, bbox_inches='tight')

