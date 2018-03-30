import os
import sys
import re

import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import config
import utils

def save(path, do_pdf=False):
    print "Writing %s" %path
    plt.savefig(path, bbox_inches='tight', dpi=300)
    if do_pdf:
        assert path.endswith('.png')
        pdf_out = re.sub('.png$', '.pdf', path)
        print "Writing %s" %pdf_out
        plt.savefig(pdf_out, bbox_inches='tight')
    return

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

