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

WHITELIST = utils.load_whitelist(config.whitelist_path)

FILTER_GTIDS = True
FILTER_GENE_PROCODE = True # only use protein coding genes

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
    cache = os.path.join(utils.get_cache_base(path), '_preproc.tsv')
    if os.path.exists(cache):
        print("Loading from %s instead of %s"%(cache, path))
        df = pd.read_csv(cache, sep='\t', index_col=0)
        psi = df.values
        strains = df.index
        gene_idx = df.columns
    else:
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

        psi, strains, gene_idx = preproc(psi, strains, gene_idx, conf_idx, ret_col_idx=ret_col_idx)
        df = pd.DataFrame(psi, index=strains, columns=gene_idx)
        if not os.path.exists(os.path.dirname(cache)): os.makedirs(os.path.dirname(cache))
        df.to_csv(cache, sep='\t')
    return psi, strains, gene_idx

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

def load_high_var_events_single_df(path, max_events, index=None, only_pc=True):
    '''Load data from each altsplice category, filter to high-var events.

        path: hdf5 path containing psi values
        map_event_to_file: dict mapping event_type -> data csv path
        max_events: filter to top max_events highest variance events
        index: filter to events in the index
        only_pc: if true, subset genes to only protein coding
    '''
    cache = utils.create_cache_base(path)
    cache = cache + '_%d_high_var_events'%max_events
    if only_pc: cache = cache + '_only_pc'
    cache = cache + '.tsv'

    if os.path.exists(cache):
        df = pd.read_csv(cache, sep='\t', index_col=0)
    else:
        print("Reading data from %s"%path)
        psi, strains, gene_idx = preproc.load_data(path)
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
        if not os.path.exists(os.path.dirname(cache)): os.makedirs(os.path.dirname(cache))
        print("Caching: %s"%cache)
        df.to_csv(cache, sep='\t')
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


