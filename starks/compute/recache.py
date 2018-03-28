import os
import sys
import h5py
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

import preproc
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import config
import utils

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
        altsplice_event_list= ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime']
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
        data = pd.read_csv(path, sep='\t', index_col=0)
        utils.cache_large_df(path.replace('.tsv', ''), data)

