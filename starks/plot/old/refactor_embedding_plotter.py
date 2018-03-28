import os
import sys

import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import utils

# input:
#    embed_dir: contains pca and tsne embeddings
#    base_plot_dir: where to write plots
#    description (e.g. altsplice-exon_skip, expression)

# tasks, fxn of embedding
#    color_cnc: color all cancer types
#    highlight_cnc: color only one cancer type (one plot per cancer)
#    tumor_normal: color tumor and normals
#    library_size: color samples by library size 
#    degscore: color samples by degragradtion score

# pca_model, pca_embeds, tsne_embeds_dict = utils.load_embeds(embed_dir)

# task: color_cnc, highlight_cnc


def _copy_dict_and_update(copy_dict, update_dict):
    ret_dict = copy_dict.copy()
    ret_dict.update(update_dict)
    return ret_dict

def cnc_plotting_args(metadata_df, highlight_cnc=None):
    '''Returns a set of objects to create figures to color cancer types.

    Input:
        metadata_df: pandas dataframe with index sample ID
                     columns for cnc and tumor status
        highlight_cnc: if None, color all cancer types,
                       otherwise color cancer types that appear in this list
                       and grey out the rest
    Returns:
        plotting_df: maps each sample ID to a color and marker type for plt.scatter
        default_kwargs: dictionary of default kwargs for plt.scatter, use if missing from plotting_df
        legend_kwargs_list: list of dictionaries specifiying each entry in the legend

    '''

    _GREY_HEX = '#D3D3D3'

    # Load color scheme
    cnc_to_color = utils.load_color_scheme(config.color_scheme_path)
    is_tumor_to_marker = {True:'o', False:'^'}

    if highlight_cnc is not None:
        # check highlight cnc are known
        if isinstance(highlight_cnc, str): highlight_cnc = [highlight_cnc]
        for cnc in highlight_cnc:
            assert cnc in cnc_to_color.keys(), "highlighted cancer type %s is unknown" %cnc
        # grey out other cancers
        for cnc, color in cnc_to_color.items():
            if cnc in highlight_cnc:
                cnc_to_color[cnc] = color
            else:
                cnc_to_color[cnc] = _GREY_HEX

    # Remove cancers that are not in the color scheme (if any)
    cancer_types = cnc_to_color.keys()
    metadata_df = metadata_df.loc[metadata_df['cnc'].isin(cancer_types)]

    # Define plotting kwargs
    default_kwargs = {'alpha': .25, 's':10}
    plotting_df = pd.DataFrame(index=metadata_df.index)
    plotting_df['facecolors'] = metadata_df['cnc'].map(cnc_to_color)
    plotting_df['marker'] = metadata_df['is_tumor'].map(is_tumor_to_marker)

    # List of kwargs used to create a legend
    legend_kwargs_list = list()
    legend_defaults = {'marker': is_tumor_to_marker[True], 'lw':0, 's': 30, 'alpha':.75, 'c':'k'}
    for cnc, color in sorted(cnc_to_color.items()):
        if highlight_cnc is not None and cnc not in highlight_cnc: continue
        legend_kwargs_list.append(_copy_dict_and_update(legend_defaults, {'label': cnc, 'c': color}))
    if highlight_cnc is not None:
        legend_kwargs_list.append(_copy_dict_and_update(legend_defaults, {'label':'Other', 'c':_GREY_HEX}))
    legend_kwargs_list.append(_copy_dict_and_update(legend_defaults, {'label': 'Tumor'}))
    legend_kwargs_list.append(_copy_dict_and_update(legend_defaults, {'label': 'Normal', 'marker': is_tumor_to_marker[False]}))

    return plotting_df, legend_kwargs_list, default_kwargs

def tumor_normal_plotting_args(metadata_df):
    '''Returns a set of objects to create figures to color tumor/normal plots

    Input:
        metadata_df: pandas dataframe with index sample ID
                     columns for cnc and tumor status
    Returns:
        plotting_df: maps each sample ID to a color and marker type for plt.scatter
        default_kwargs: dictionary of default kwargs for plt.scatter, use if missing from plotting_df
        legend_kwargs_list: list of dictionaries specifiying each entry in the legend

    '''

    is_tumor_to_color = {True: 'lightblue', False: 'orange'}
    default_kwargs = {'alpha': .25, 's':10}
    is_tumor_to_marker = {True:'o', False:'^'}

    plotting_df = pd.DataFrame(index=metadata_df.index)
    plotting_df['facecolors'] = metadata_df['is_tumor'].map(is_tumor_to_color)
    plotting_df['marker'] = metadata_df['is_tumor'].map(is_tumor_to_marker)

    # List of kwargs used to create a legend
    legend_kwargs_list = list()
    legend_defaults = {'lw':0, 's': 30, 'alpha':.75}
    tumor_args  = {'label': 'Tumor',  'c': is_tumor_to_color[True],  'marker': is_tumor_to_marker[True]}
    normal_args = {'label': 'Normal', 'c': is_tumor_to_color[False], 'marker': is_tumor_to_marker[False]}
    legend_kwargs_list.append(_copy_dict_and_update(legend_defaults, tumor_args))
    legend_kwargs_list.append(_copy_dict_and_update(legend_defaults, normal_args))
    return plotting_df, legend_kwargs_list, default_kwargs

def _unpack_scatter_kwargs(marker_kwargs, plotting_df_subset):
    '''Unpack cols from plotting_df_subset into marker_kwargs dict.
    '''
    for kw in plotting_df_subset.columns:
        if kw == 'marker': continue
        marker_kwargs[kw] = plotting_df_subset[kw].values
    return

def plot_embeddings(embed_df, plotting_df, default_kwargs=None, legend_kwargs_list=None, title=None, ax=None):
    '''Create scatter plot of embeddings

    Input:
        embed_df: samples x embed_dim dataframe
        plotting_df: samples x kwargs (requires 'marker' column)
        default_kwargs: dict of plt.scatter kwargs, uses these values if missing from plotting_df
        legend_kwargs_list: list of dict, each entry defines an entry in the legend
        title: title for the plot
        ax: plt object to house plot

    Returns
        fig, ax: matplotlib objects for plot

    '''

    if ax is None: fig, ax = plt.subplots()
    else: fig = ax.get_figure()
    if default_kwargs is None: default_kwargs = dict()

    # Scatter requires one call for each marker type
    assert 'marker' in plotting_df
    marker_groups = plotting_df.groupby(plotting_df['marker'])
    for marker, plotting_df_subset in marker_groups:
        embeds = embed_df.loc[plotting_df_subset.index].values
        scatter_kwargs = default_kwargs.copy()
        _unpack_scatter_kwargs(scatter_kwargs, plotting_df_subset)
        scatter_kwargs['marker'] = marker
        ax.set_xticks([]); ax.set_yticks([])
        ax.scatter(embeds[:, 0], embeds[:, 1], **scatter_kwargs)

    if legend_kwargs_list is not None:
        add_legend(ax, legend_kwargs_list, loc='bottom', ncol=5)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    if title is not None: ax.set_title(title)
    return fig, ax

def plot_continuous_color_embeddings(embed_df, scale_series, axis_title=None, cbar_title=None, ax=None):
    if ax is None: fig, ax = plt.subplots()
    else: fig = ax.get_figure()
    indx = embed_df.index.intersection(scale_series.index)
    scatter = ax.scatter(embed_df.loc[indx].iloc[:, 0], embed_df.loc[indx].iloc[:, 1], c=scale_series.loc[indx].values, cmap='Blues')
    cbar = fig.colorbar(scatter)

    if axis_title is not None: ax.set_title(axis_title)
    if cbar_title is not None:
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel(cbar_title, rotation=270)
    return fig, ax

def add_legend(ax, label_kwargs_list, loc='right', ncol=1):
    '''Adds legend to ax

    ax: axis object to anchor label
    label_kwargs_list: list of plotting kwarg dicts
    '''
    labels = list()
    scatters = list()
    for kwargs in label_kwargs_list:
        assert 'label' in kwargs
        labels.append(kwargs['label'])
        scatters.append(ax.scatter([],[], **kwargs))

    if loc == 'right':
        bbox_to_anchor = (1.05, 1.0)
        loc = 'upper left'
    elif loc == 'bottom':
        bbox_to_anchor = (.5, -.1)
        loc = 'upper center'
    else:
        raise ValueError

    ax.legend(scatters, labels,
              bbox_to_anchor=bbox_to_anchor,
              loc=loc,
              borderaxespad=0.0,
              ncol=ncol)
    return

