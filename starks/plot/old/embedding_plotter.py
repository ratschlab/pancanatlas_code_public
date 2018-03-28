import os
import sys

import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
import utils


def load_subtype_color_tumor_marker_kwargs(metadata_df, subtype_name, highlight_cnc):
    _DEFAULT_COLOR = '#D3D3D3' # grey
    _UNK_COLOR = '#ffc0cb' # peach pink

    default_kwargs = {'alpha': .25, 's':10}
    plotting_df = pd.DataFrame(index=metadata_df.index)
    plotting_df['facecolors'] = _DEFAULT_COLOR
    plotting_df['facecolors'].loc[metadata_df['cnc'] == highlight_cnc] = _UNK_COLOR
    is_tumor_to_marker = {True:'o', False:'^'}
    plotting_df['marker'] = metadata_df['is_tumor'].map(is_tumor_to_marker)

    # Add the subtype information
    subtype_series = metadata_df[subtype_name].loc[metadata_df['cnc'] == highlight_cnc]
    subtype_series = subtype_series.dropna()
    subtype_series = subtype_series.loc[subtype_series.map(lambda x: not 'normal' in x.lower())]
    if subtype_series.size == 0: return None, None, None
    values = np.sort(np.append(subtype_series.unique(), 'Normal'))
    colors = sns.color_palette("Set1", n_colors=values.size, desat=.5).as_hex()
    color_lut = dict(zip(values, colors))
    plotting_df['facecolors'].loc[subtype_series.index] = subtype_series.map(lambda x: color_lut[x])
    # add colors for normal samples
    plotting_df['facecolors'].loc[np.logical_and(metadata_df['cnc'] == highlight_cnc, ~metadata_df['is_tumor'])] = color_lut['Normal']

    # List of kwargs used to create a legend
    legend_kwargs_list = list()
    legend_defaults = {'marker': 'o', 'lw':0, 's': 30, 'alpha':.75}
    for stype, color in sorted(color_lut.items()):
        if stype == 'Normal':
            update_dict = {'c':color, 'label': stype, 'marker': '^'}
        else:
            update_dict = {'c':color, 'label': stype}
        legend_kwargs_list.append(_copy_and_update(legend_defaults, update_dict))
    legend_kwargs_list.append(_copy_and_update(legend_defaults, {'c':_DEFAULT_COLOR, 'label':'Other'}))
    if _UNK_COLOR in plotting_df['facecolors'].values:
        legend_kwargs_list.append(_copy_and_update(legend_defaults, {'c':_UNK_COLOR, 'label':'Subtype Unknown'}))
    return plotting_df, legend_kwargs_list, default_kwargs

def load_cnc_color_tumor_marker_kwargs(metadata_df, highlight_cnc=None):
    # Load color scheme
    cnc_to_color = utils.load_color_scheme(config.color_scheme_path)
    is_tumor_to_marker = {True:'o', False:'^'}

    if highlight_cnc is not None:
        # check highlight cnc's
        if isinstance(highlight_cnc, str): highlight_cnc = [highlight_cnc]
        for cnc in highlight_cnc:
            assert cnc in cnc_to_color.keys(), "highlighted cancer type %s is unknown" %cnc

        # grey out colors
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
    legend_defaults = {'marker': is_tumor_to_marker[True], 'lw':0, 's': 30, 'alpha':.75}
    for cnc, color in sorted(cnc_to_color.items()):
        if highlight_cnc is not None and cnc not in highlight_cnc: continue
        legend_kwargs_list.append(_copy_and_update(legend_defaults, {'c':color, 'label': cnc}))
    if highlight_cnc is not None:
        legend_kwargs_list.append(_copy_and_update(legend_defaults, {'c':_GREY_HEX, 'label':'Other'}))
    legend_kwargs_list.append(_copy_and_update(legend_defaults, {'c':'k', 'label': 'Tumor'}))
    legend_kwargs_list.append(_copy_and_update(legend_defaults,
        {'c':'k', 'label': 'Normal', 'marker': is_tumor_to_marker[False]}))

    return plotting_df, legend_kwargs_list, default_kwargs

def process(plot_dir, pca_model, pca_embeds, tsne_embeds_dict,
            plotting_df, legend_kwargs_list, default_kwargs,
            do_pca=True, do_tsne=True,
            tsne_pp_plot_set='all', tsne_lr_plot_set='all'):

    '''Default plotting routine for pca/tsne data.

    Plots 3x3 pca scatter plot to pca_embeds.png
    Plots tsne plots to tsne_embeds_pp_%pp_lr_%lr.png

    pca_model: instance of sklearn pca model (for explained variance)
    pca_embeds: data frame of pca embeddings
    tsne_embeds_dict: dict with keys (pp,lr) values tsne embeddings
    plotting_df: dataframe with ax.scatter kwargs
    plot_dir: where to save plots
    do_pca: if true, write pca plot
    do_tsne: if true, write tsne plots
    tsne_pp_plot_set: iterable, subset to only tsne runs with given pp values
    tsne_lr_plot_set: iterable, subset to only tsne runs with given lr values
    '''

    if do_pca:
        pca_plot_out = os.path.join(plot_dir, 'pca_embeds.png')
        fig, axes = plot_pca_embeds(pca_embeds, plotting_df, **default_kwargs)
        add_legend(axes[-1][1], legend_kwargs_list, loc='bottom', ncol=5)
        add_pca_variance_expl(axes, pca_model)
        fig.suptitle('Raw Expression PCA')
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        print "Writing %s" %pca_plot_out
        plt.savefig(pca_plot_out, bbox_inches='tight', dpi=300)
        plt.close()

    if do_tsne:
        for (pp, lr), tsne_embeds in tsne_embeds_dict.items():
            if not tsne_pp_plot_set == 'all' and not pp in tsne_pp_plot_set: continue
            if not tsne_lr_plot_set == 'all' and not lr in tsne_lr_plot_set: continue
            tsne_plot_out = os.path.join(plot_dir, 'tsne_embeds_pp_%d_lr_%d.png' %(pp, lr))

            fig, ax = plot_tsne_embeds(tsne_embeds_dict[(pp,lr)], plotting_df, **default_kwargs)
            add_legend(ax, legend_kwargs_list, loc='bottom', ncol=5)
            fig.suptitle('Raw Expression tSNE(perplexity=%d, learning_rate=%d)' %(pp, lr))
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            print "Writing %s" %tsne_plot_out
            plt.savefig(tsne_plot_out, bbox_inches='tight', dpi=300)
            plt.close()

def _copy_and_update(copy_dict, update_dict):
    ret_dict = copy_dict.copy()
    ret_dict.update(update_dict)
    return ret_dict

def add_legend(ax, label_kwargs_list, loc='right', ncol=1):
    '''Adds legend to ax

    ax: axis object to anchor label
    label_kwargs_list: list of plotting kwarg dicts
    cmap: used to extract color from kwarg dicts

    if loc == right, give the top right axis
    if loc == bottom, give the bottom center axis
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

def auto_format_axes(naxes):
    if naxes == 1:
        fig, axes = plt.subplots()
    elif naxes == 2:
        fig, axes = plt.subplots(1,2)
    elif naxes <= 4:
        fig, axes = plt.subplots(2,2)
    elif naxes <= 9:
        fig, axes = plt.subplots(3,3)
    else:
        raise NotImplementedError
    return fig, axes

def _update_marker_kwargs(marker_kwargs, plotting_df_subset):
    '''Unpack cols from plotting_df_subset into marker_kwargs dict.
    '''
    for kw in plotting_df_subset.columns:
        if kw == 'marker': continue
        marker_kwargs[kw] = plotting_df_subset[kw].values
    return

def _scatter_all_axes(axes, data, **kwargs):
    '''Scatter slices of data for each dim in axes
    '''
    for i in range(axes.shape[0]):
        for j in range(axes.shape[1]):
            ax = axes[i][j]
            ax.set_xticks([])
            ax.set_yticks([])
            if i == j: continue
            scat = ax.scatter(data[:, i], data[:, j], **kwargs)
    return

def add_pca_embeddings(axes, embeds, **kwargs):
    '''Adds embedding slices to each axis.
    '''
    if isinstance(axes, plt.Axes):
        scat = axes.scatter(embeds[:, 0], embeds[:, 1], **kwargs)
    else:
        for i in range(axes.shape[0]):
            for j in range(axes.shape[1]):
                ax = axes[i][j]
                ax.set_xticks([])
                ax.set_yticks([])
                if i == j: continue
                scat = ax.scatter(embeds[:, i], embeds[:, j], **kwargs)
    return scat

def add_pca_variance_expl(axes, model):
    '''Adds a text box to the diagonal with percent variance explained.
    '''
    for i in range(axes.shape[0]):
        for j in range(axes.shape[1]):
            if i != j: continue
            ratio = model.explained_variance_ratio_[i]
            ax = axes[i][j]
            ax.text(.5, .5, "PC %d\n(%.2f%%)"%(i+1, 100*ratio),
                            horizontalalignment='center', verticalalignment='center',
                            transform=ax.transAxes, fontsize=10)
    return

def plot_pca_embeds(embed_df, plotting_df, axes=None, **default_kwargs):
    '''Plot first 3 PCs against each other

    embed_df: dataframe of embeddings
    plotting_df: dataframe of scatter kwargs specific to each data point
    kwargs: arguments for plt.scatter
    '''
    if axes is None:
        fig, axes = plt.subplots(3,3)
    else:
        assert axes.ravel().size == 9
        fig = axes.ravel()[0].get_figure()

    assert 'marker' in plotting_df
    marker_groups = plotting_df.groupby(plotting_df['marker'])
    for marker, plotting_df_subset in marker_groups:
        embeds = embed_df.loc[plotting_df_subset.index].values
        scatter_kwargs = default_kwargs.copy()
        _update_marker_kwargs(scatter_kwargs, plotting_df_subset)
        scatter_kwargs['marker'] = marker
        _scatter_all_axes(axes, embeds, **scatter_kwargs)
    return fig, axes


def plot_tsne_embeds(embed_df, plotting_df, ax=None, **kwargs):
    '''Plots tsne embeds in a single axis.

    embed_df: dataframe of embeddings
    plotting_df: dataframe of scatter kwargs specific to each data point
    kwargs: arguments for plt.scatter
    '''
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()
    assert 'marker' in plotting_df
    marker_groups = plotting_df.groupby(plotting_df['marker'])
    for marker, plotting_df_subset in marker_groups:
        embeds = embed_df.loc[plotting_df_subset.index].values
        scatter_kwargs = kwargs.copy()
        _update_marker_kwargs(scatter_kwargs, plotting_df_subset)
        scatter_kwargs['marker'] = marker
        ax.set_xticks([])
        ax.set_yticks([])
        ax.scatter(embeds[:, 0], embeds[:, 1], **scatter_kwargs)
    return fig, ax

def get_embeddings(embed_dir):
    pca_model, pca_embeds, tsne_embeds_dict = utils.load_embeds(embed_dir)
    metadata_df = utils.load_metadata_df(config.metadata_path, pca_embeds.index)
