import os
import sys
import argparse

import numpy as np
import pandas as pd

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc

sys.path.append('../..')
import config
import utils
import compute.alt_splice as preproc

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

def load_high_var_events_single_df(path, cache_base, max_events=None, index=None, only_pc=False):
    '''Load data from each altsplice category, filter to high-var events.

        path: tsv input for embeddings
        cache_base: base for loading/saving data
        map_event_to_file: dict mapping event_type -> data csv path
        max_events: filter to top max_events highest variance events
        index: filter to events in the index
        only_pc: if true, subset genes to only protein coding
    '''
    if only_pc: cache_path = cache_base + '_only_pc'
    cache_path = cache_path + '_%d_high_var_events.tsv'%max_events
    if os.path.exists(cache_path):
        df = pd.read_csv(cache_path, sep='\t', index_col=0)

    print("Reading data from %s"%path)
    df = pd.read_csv(path)
    psi = df.values
    strains = df.index
    gene_idx = df.columns
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

def _filter_to_protein_coding(psi, gene_idx, gene_names):
    pc_ensgs = np.loadtxt(config.fn_protein_coding_list, delimiter='\t', dtype=str, usecols=[0])
    import ipdb; ipdb.set_trace()
    return psi, gene_idx

def filter_to_high_var(data, columns, nkeep):
    '''Filter to the top nkeep high variance columns
    '''
    if nkeep is None: return data, columns
    if nkeep <=1: nkeep = int(data.shape[1] * nkeep)
    var = np.var(data, axis=0)
    assert var.size == data.shape[1]
    keep_cols = np.argsort(var)[-nkeep:]
    data = data[:, keep_cols]
    columns = columns[keep_cols]
    return data, columns


