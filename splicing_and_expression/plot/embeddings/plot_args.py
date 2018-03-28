import os
import sys
import pandas as pd
import numpy as np

import matplotlib; matplotlib.use('Agg')
import seaborn as sns


sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import config
import utils

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

def load_subtype_color_tumor_marker_kwargs(metadata_df, subtype_name, highlight_cnc):
    _DEFAULT_COLOR = '#D3D3D3' # grey
    _UNK_COLOR     = '#ffc0cb' # peach pink

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
        legend_kwargs_list.append(_copy_dict_and_update(legend_defaults, update_dict))
    legend_kwargs_list.append(_copy_dict_and_update(legend_defaults, {'c':_DEFAULT_COLOR, 'label':'Other'}))
    if _UNK_COLOR in plotting_df['facecolors'].values:
        legend_kwargs_list.append(_copy_dict_and_update(legend_defaults, {'c':_UNK_COLOR, 'label':'Subtype Unknown'}))

    return plotting_df, legend_kwargs_list, default_kwargs


