import os
import sys
import h5py
import numpy as np
import scipy as sp
import pandas as pd

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

def preproc(psi, strains, gene_idx, conf_idx):
    '''Mask gene_idx to those in conf_idx.
    '''
    assert psi.shape == (strains.size, gene_idx.size)
    gene_idx = gene_idx[conf_idx]
    psi = psi[:, conf_idx]
    psi, row_mask, col_mask = replace_nans(psi)
    col_idx = conf_idx[col_mask]
    strains = strains[row_mask]
    gene_idx = gene_idx[col_mask]
    assert psi.shape == (strains.size, gene_idx.size)
    return psi, strains, gene_idx, col_idx


def load_data(path):
    data = h5py.File(path, 'r')
    psi = data['psi'][:]
    print "\t psi.ngbytes = %.1f" %(psi.nbytes * 1e-9)
    gene_idx = data['gene_idx'][:]
    conf_idx = data['conf_idx'][:]
    strains = data['strains'][:]
    data.close()
    print

    psi, strains, gene_idx, col_idx = preproc(psi, strains, gene_idx, conf_idx)
    return psi, strains, gene_idx, col_idx

def clean_strain(strain):
    strain = np.array([x.replace('.aligned', '').replace('.npz', '') for x in strain])
    return strain


