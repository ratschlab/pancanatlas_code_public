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

if __name__ == '__main__':
    do_pca = True
    do_tsne = True
    tsne_pp_plot_set = 'all' #[100]
    tsne_lr_plot_set = 'all' #[500]

    # Load embeds, metadata
    embeds_dir = os.path.join(config.embed_dir, 'raw_expression')
    pca_model, pca_embeds, tsne_embeds_dict = utils.load_embeds(embeds_dir)
    metadata_df = utils.load_metadata_df(config.metadata_path, pca_embeds.index)

    plot_dir = os.path.join(config.plot_dir, 'raw_expression')
    if not os.path.exists(plot_dir): os.makedirs(plot_dir)
    embedding_plotter.process(pca_model, pca_embeds, tsne_embeds_dict,
                              metadata_df, plot_dir,
                              do_pca=True do_tsne=True,
                              tsne_pp_plot_set=tsne_pp_plot_set,
                              tsne_lr_plot_set=tsne_lr_plot_set)


