import os
import sys
import numpy as np
import scipy as sp

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import config

def print_nan_stats(data, data_is_fin=None):
    '''Prints a table showing how many rows/cols left at each thold.
    '''
    if data_is_fin is None: data_is_fin = np.isfinite(data)
    rows_pct_fin = np.mean(data_is_fin, 1)
    assert rows_pct_fin.size == data_is_fin.shape[0]
    cols_pct_fin = np.mean(data_is_fin, 0)
    assert cols_pct_fin.size == data_is_fin.shape[1]
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

def impute_non_finite(data, data_is_fin=None):
    '''Replace non-finite entries with mean of columns.
    '''
    if data_is_fin is None: data_is_fin = np.isfinite(data)
    col_means = sp.nanmean(data, 0)
    nan_ridx, nan_cidx = np.where(np.invert(data_is_fin))
    data[nan_ridx, nan_cidx] = col_means[nan_cidx]

def replace_nans(data, col_thold=.7, row_thold=.9):
    '''Replace nan entries in each row with mean of column
    '''
    data_is_fin = np.isfinite(data)
    print "Stats before filtering"
    print_nan_stats(data, data_is_fin)

    print "Removing columns with < %d%% finite values" %int(100 * col_thold)
    cols_pct_fin = np.mean(data_is_fin, 0)
    assert cols_pct_fin.size == data_is_fin.shape[1]
    col_mask = cols_pct_fin > col_thold
    data = data[:, col_mask]
    data_is_fin = data_is_fin[:, col_mask]
    print_nan_stats(data, data_is_fin)
    print ""

    print "Removing rows with < %d%% finite values" %int(100 * row_thold)
    rows_pct_fin = np.mean(data_is_fin, 1)
    row_mask = rows_pct_fin > row_thold
    data = data[row_mask]
    data_is_fin = data_is_fin[row_mask]
    print_nan_stats(data, data_is_fin)
    print ""

    print "Replacing nans with average of column"
    col_means = sp.nanmean(data, 0)
    nan_ridx, nan_cidx = np.where(np.invert(data_is_fin))
    data[nan_ridx, nan_cidx] = col_means[nan_cidx]
    return data, row_mask, col_mask

def filter_sample_ids(sample_ids, counts=None):
    #TODO: filter samples that have unusually low gene counts
    sample_ids_mask = sample_ids != 'NA'
    print "Missing %d sample_ids" %np.invert(sample_ids_mask).sum()
    return sample_ids_mask


def filter_genes(gene_ids, counts=None, only_protein_coding=True, avg_lbound=5):
    '''Returns a mask over gene_ids according to filter conditions.

    gene_ids: numpy array of gene names (i.g. ENSG...)
    counts: genes x samples
    only_protein_coding (bool): filter out non-protein coding genes
    avg_lbound (int): filter out genes with avg count < avg_lbound
    '''
    gene_mask = np.ones(gene_ids.size, dtype=bool)
    if only_protein_coding:
        protein_coding = config.load_protein_coding_genes()
        pc_mask = np.in1d(gene_ids, protein_coding)
        gene_mask = np.logical_and(gene_mask, pc_mask)
        print "Filtering to %d protein coding genes" %pc_mask.sum()

    if avg_lbound > 0:
        assert not counts is None
        count_mask = counts.mean(1) > avg_lbound
        gene_mask = np.logical_and(gene_mask, count_mask)
        print "Filtering to %d genes with avg counts_raw > %d" %(gene_mask.sum(), avg_lbound)

    return gene_mask

