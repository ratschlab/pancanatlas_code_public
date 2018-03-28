import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

_ALPHA = 1

def _process_colors(kwargs_config, cmap):
    if 'color_index' in kwargs_config:
        assert 'ndata' in kwargs_config
        color_index = kwargs_config.pop('color_index')
        ndata = kwargs_config.pop('ndata')
        ncats = kwargs_config.pop('ncats', None)
        if not ncats is None:
            color_index = float(color_index) / ncats
        else:
            color_index = color_index % cmap.N
        return np.tile(cmap(color_index), (ndata, 1))

def unpack_kwargs_hollow(kwargs_config, cmap):
    kwargs_config = kwargs_config.copy()
    kwargs = {'facecolors':'none',
              'alpha': _ALPHA,
              's': 20,
              'lw': 1}
    if isinstance(cmap, basestring): cmap = plt.get_cmap(cmap)
    color_tile = _process_colors(kwargs_config, cmap)
    kwargs_config['edgecolors'] = color_tile
    kwargs.update(kwargs_config)
    return kwargs

def unpack_kwargs_meta(kwargs_config, cmap):
    kwargs_config = kwargs_config.copy()
    if isinstance(cmap, basestring): cmap = plt.get_cmap(cmap)
    kwargs_config.pop('is_meta', None)
    kwargs = {'alpha': _ALPHA,
              'facecolors': 'black',
              's': 30,
              'lw': 1}

    color_tile = _process_colors(kwargs_config, cmap)
    kwargs_config['edgecolors'] = color_tile
    kwargs.update(kwargs_config)
    return kwargs

def unpack_kwargs_filled(kwargs_config, cmap):
    kwargs_config = kwargs_config.copy()
    if isinstance(cmap, basestring): cmap = plt.get_cmap(cmap)
    kwargs = {'alpha': _ALPHA,
#              'edgecolors': 'black',
#              'lw': 1,
              's': 20}

    color_tile = _process_colors(kwargs_config, cmap)
    kwargs_config['facecolors'] = color_tile
    kwargs.update(kwargs_config)
    return kwargs

def unpack_kwargs_legend(kwargs, cmap):
    kwargs = kwargs.copy()
    if 'ndata' in kwargs: kwargs['ndata'] = 1
    kwargs = unpack_kwargs_filled(kwargs, cmap)
    kwargs['lw'] = 0
    kwargs['s'] = 20
    kwargs['alpha'] = .75
    return kwargs


def get_cancer_type_only_maps(cnc_series, default_kwargs, max_nlabels=40):
    map_label_to_index = dict()
    map_label_to_kwargs = dict()
    legend_kwargs_list = list()

    groups = sorted(cnc_series.groupby(cnc_series), key=lambda x: -x[1].size)[:max_nlabels]

    for color_index, (label, data) in enumerate(sorted(groups)):
        index = data.index

        legend_kwargs = dict()
        legend_kwargs['label'] = "%s [%d]" %(label, index.size)
        legend_kwargs['color_index'] = color_index
        legend_kwargs['ndata'] = 1
        legend_kwargs_list.append(legend_kwargs)

        map_label_to_index[label] = index
        map_label_to_kwargs[label] = default_kwargs(color_index, index.size)

    return map_label_to_index, map_label_to_kwargs, legend_kwargs_list


def get_label_submarker_maps(cnc_series, is_tumor_series, default_kwargs, max_nlabels=40):
    '''Configures plotting arguments, color for cancer type, marker for tumor/normal.

    cnc_series: used to choose color
    is_tumor_series: used to choose marker
    default_kwargs: function of (color_index, size), sets other default kwargs.
    max_nlabels: only keep the top most common labels

    One label for each color (cancer type)
    Two meta labels in black for tumor/normal marker
    '''
    tumor_marker = 'o'
    normal_marker = '^'
    map_label_to_index = dict()
    map_label_to_kwargs = dict()
    legend_kwargs_list = list()

    groups = sorted(cnc_series.groupby(cnc_series), key=lambda x: -x[1].size)[:max_nlabels]

    for color_index, (label, data) in enumerate(sorted(groups)):
        index = data.index
        is_tumor = is_tumor_series.loc[index].values
        tumor_label = label + 'Tumor'
        normal_label = label + 'Normal'

        legend_kwargs = dict()
        legend_kwargs['label'] = "%s [%d]" %(label, index.size)
        legend_kwargs['color_index'] = color_index
        legend_kwargs['ndata'] = 1
        legend_kwargs_list.append(legend_kwargs)

        if is_tumor.sum() > 0:
            map_label_to_index[tumor_label] = index[is_tumor]
            map_label_to_kwargs[tumor_label] = default_kwargs(color_index, is_tumor.sum())
            map_label_to_kwargs[tumor_label]['marker'] = tumor_marker

        is_normal = np.invert(is_tumor)
        if is_normal.sum() > 0:
            map_label_to_index[normal_label] = index[is_normal]
            map_label_to_kwargs[normal_label] = default_kwargs(color_index, is_normal.sum())
            map_label_to_kwargs[normal_label]['marker'] = normal_marker

    legend_kwargs_list.append({'label': 'Tumor',
                              'marker': tumor_marker,
                              'c': 'k'})

    legend_kwargs_list.append({'label': 'Normal',
                              'marker': normal_marker,
                              'c': 'k'})

    return map_label_to_index, map_label_to_kwargs, legend_kwargs_list

def test_label_maps(metadata_df, map_keys_to_kwargs, default_kwargs):
    '''
    '''
