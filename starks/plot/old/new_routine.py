import os
import sys
import numpy as np
import pandas as pd
import h5py

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import config
import utils

if __name__ == '__main__':
    # Load embeds, metadata
    embeds_dir = os.path.join(config.embed_dir, 'raw_expression')
    pca_model, pca_embeds, tsne_embeds_dict = utils.load_embeds(embeds_dir)
    metadata_df = utils.load_metadata_df(config.metadata_path, pca_embeds.index)

    # Remove rare cancer types (just for color scheme)
    max_cancer_types = 20
    cancer_types = np.sort(metadata_df['cnc'].value_counts().index[:max_cancer_types])
    metadata_df = metadata_df.loc[metadata_df['cnc'].isin(cancer_types)]
    joint_index = metadata_df.index & pca_embeds.index
    metadata_df = metadata_df.loc[joint_index]

    # Define a color scheme
    cmap = plt.get_cmap('tab20')
    color_list = [cmap(i) for i in range(20)[::2]] + [cmap(i) for i in range(20)[1::2]]
    cnc_to_color = dict(zip(cancer_types, color_list))
    is_tumor_to_marker = {True:'o', False:'^'}

    # Define plotting kwargs
    default_kwargs = {'alpha': .25, 's':2}
    plotting_df = pd.DataFrame(index=joint_index)
    plotting_df['c'] = metadata_df['cnc'].map(cnc_to_color)
    plotting_df['marker'] = metadata_df['is_tumor'].map(is_tumor_to_marker)

    # Legend kwargs, label each cancer type and tumor/normal shapes
    legend_kwarg_list = list()
    for cnc, color in cnc_to_color.items():
        legend_kwargs_list.append({'c':color, 'label': cnc, 'marker': is_tumor_marker[True]})
    legend_kwargs_list.append({'label': 'Tumor','marker': is_tumor_to_marker[True], 'c': 'k'})
    legend_kwargs_list.append({'label': 'Normal', 'marker': is_tumor_to_marker[False], 'c': 'k'})

    pca_fig, pca_axes = plt.subplots(3,3)
    pca_plot(pca_embeds, plotting_df, legend_kwargs_list, **default_kwargs)
