import os
import sys
import h5py
import numpy as np
import scipy as sp
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import config
import utils
import preproc
import embed


DO_PCA = True
DO_TSNE = True and DO_PCA

PCA_NCOMPONENTS = 100
TSNE_LR_LIST = [100, 500, 1000]
TSNE_PP_LIST = [20, 40, 50, 75, 100, 200, 500]

#Common variants: > 3% (non-2 values)
_DEBUG = True

def stream_filter_uncommon_variants(gt_handle, common_freq=.97, missing_freq=.05):
    '''Read in data by cols, remove columns with too many missing reference values.

    rows are genes, cols are samples
    streams data by row
    remove rows that have:
        percentage of samples (reference + missing) is > common_freq
        percentage of samples missing is > missing_freq
    '''
    mask = list()
    data = list()

    assert gt_handle.shape[0] > gt_handle.shape[1] # genes x samples
    full_gb = np.prod(gt_handle.shape) * 8 * 1e-9
    chunk_size = 100 * gt_handle.chunks[0]
    indices = np.arange(0, gt_handle.shape[0], step=chunk_size)
    for index, start in enumerate(indices):
        if index == indices.size - 1:
            chunk_data = gt_handle[start:, :]
        else:
            chunk_data = gt_handle[start:indices[index+1], :]

        # Keep gene if samples meet:
        #   less than 97% are reference (filter rare alterations)
        #   less than 5% of samples are nan (missing data)
        pct_reference = np.mean(chunk_data == 2, axis=1)
        pct_nan = np.mean(np.isnan(chunk_data), axis=1)
        mask_to_commonly_altered_genes = (pct_reference + pct_nan) < common_freq
        mask_to_genes_with_few_missing = pct_nan < missing_freq
        chunk_mask = mask_to_commonly_altered_genes * mask_to_genes_with_few_missing

        mask.extend(chunk_mask.tolist())
        if len(chunk_mask) > 0 and np.any(chunk_mask):
            chunk_data = chunk_data[chunk_mask]
            data.extend(chunk_data.tolist())
        mask_mean = np.mean(mask) if len(mask) > 0 else 0
        filler = (index, indices.size, mask_mean, full_gb*mask_mean)
        print "%10d of %10d Mean mask: %.2f, Est Gb: %.1f\r"%filler,

    print ""
    print "Filtered genes if \n  > %d%% of values were (nan OR ref)\n  > %d%% values were nan"%(100*common_freq, 100*missing_freq)
    print "%d -> %d genes" %(len(mask), sum(mask))
    data = np.array(data)
    mask = np.array(mask)
    preproc.impute_non_finite(chunk_data)
    return data, mask

def load_data(config):
    '''Load the germline & somatic dataframes.
    '''

    germline_data = h5py.File(config.genotype_germline_path, 'r')
    index = germline_data['gtid'][:]
    germline_gt, germline_pos_mask = stream_filter_uncommon_variants(germline_data.get('gt'))
    if germline_gt.shape[0] != index.size: germline_gt = germline_gt.T
    germline_pos = germline_data['pos'][:][germline_pos_mask]
    germline_pos = np.array(['_'.join(x.astype(str)) for x in germline_pos])
    assert germline_gt.shape == (index.size, germline_pos.size)
    germline_data.close()
    germline_df = pd.DataFrame(germline_gt, index=index, columns=germline_pos)

    somatic_data = h5py.File(config.genotype_somatic_path, 'r')
    somatic_gt = somatic_data['gt'][:]
    if somatic_gt.shape[0] != index.size: somatic_gt = somatic_gt.T
    somatic_pos = somatic_data['pos'][:]
    assert somatic_gt.shape == (index.size, somatic_pos.size)
    somatic_data.close()
    somatic_df = pd.DataFrame(somatic_gt, index=index, columns=somatic_pos)
    return somatic_df, germline_df

CACHE = True
CACHE_DIR = os.path.expanduser('~/cache')
if not os.path.exists(CACHE_DIR): os.makedirs(CACHE_DIR)
SOMATIC_CACHE = os.path.join(CACHE_DIR, 'somatic_genotype_df.tsv')
GERMLINE_CACHE = os.path.join(CACHE_DIR, 'germline_genotype_df.tsv')

if __name__ == '__main__':
    if CACHE and os.path.exists(SOMATIC_CACHE) and os.path.exists(GERMLINE_CACHE):
        print "Reading from cache"
        somatic_df = pd.read_csv(SOMATIC_CACHE, sep='\t', index_col=0)
        germline_df = pd.read_csv(GERMLINE_CACHE, sep='\t', index_col=0)
    else:
        somatic_df, germline_df = load_data(config)
        if CACHE:
            somatic_df.to_csv(SOMATIC_CACHE, sep='\t')
            germline_df.to_csv(GERMLINE_CACHE, sep='\t')

    for name, df in zip(['somatic', 'germline'], [somatic_df, germline_df]):
        embed_dir = os.path.join(config.embed_dir,'genotype', name)
        if not os.path.exists(embed_dir): os.makedirs(embed_dir)
        embed.run_embed_proc(df.values, df.index, embed_dir,
                             pca_ncomponents=PCA_NCOMPONENTS,
                             tsne_lr_list=TSNE_LR_LIST,
                             tsne_pp_list=TSNE_PP_LIST,
                             do_tsne=DO_TSNE)
