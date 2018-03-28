import os
import sys
import numpy as np
import pandas as pd

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import embedding_plotter

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import config
import utils

#TODO: something weird with the metadata/embedding indices here..
if __name__ == '__main__':
    do_pca = True
    do_tsne = True
    tsne_pp_plot_set = [100]
    tsne_lr_plot_set = [500]

    plot_base_dir = os.path.join(config.plot_dir, 'genotype')
    embed_base_dir = os.path.join(config.embed_dir, 'genotype')

    gtype_list = ['germline', 'somatic']
    for gtype in gtype_list:
        embed_dir = os.path.join(embed_base_dir, gtype)
        assert os.path.exists(embed_dir)
        plot_dir = os.path.join(plot_base_dir, gtype)
        if not os.path.exists(plot_dir): os.makedirs(plot_dir)

        pca_model, pca_embeds, tsne_embeds_dict = utils.load_embeds(embed_dir)
        metadata_df = utils.load_metadata_df(config.metadata_path)
        embedding_plotter.process(pca_model, pca_embeds, tsne_embeds_dict,
                                  metadata_df, plot_dir,
                                  do_pca=do_pca, do_tsne=do_tsne,
                                  tsne_pp_plot_set=tsne_pp_plot_set,
                                  tsne_lr_plot_set=tsne_lr_plot_set)
