import os
import sys
import re
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
from matplotlib import colors
import seaborn as sns
import matplotlib.pyplot as plt

def test_data():
    nsamples = 1000
    nfeatures = 1500
    tumor = np.random.randn(nsamples, nfeatures) * .5 - .5
    normal = np.random.randn(nsamples, nfeatures) * .5 + .5

    row_labels = ['tumor'] * nsamples + ['normal'] * nsamples
    base_color = sns.color_palette('Set2', 1)[0]
    tumor_color = base_color + (1,)
    normal_color = base_color + (.5,)
    map_label_to_color = {'tumor': tumor_color, 'normal': normal_color}
    row_colors = map(lambda x: map_label_to_color[x], row_labels)

    data = np.vstack((tumor, normal))
    assert data.shape == (2 * nsamples, nfeatures)
    return data, row_colors

def _override_sns_row_colors(graph, row_colors):

    if isinstance(row_colors[0], tuple):
        # row_colors are in rgb(a) form
        try:
            row_colors = row_colors.tolist()
        except AttributeError:
            pass
        unq_colors, color_class = np.unique(row_colors, axis=0, return_inverse=True)
        unq_colors = map(lambda x: tuple(x), unq_colors)
    else:
        unq_colors, color_class = np.unique(row_colors, return_inverse=True)
        unq_colors = unq_colors.tolist()

    rcax = graph.ax_row_colors
    rcax.clear()
    cmap = colors.ListedColormap(unq_colors)
    rcax.imshow(np.matrix(color_class).T, aspect='auto', cmap=cmap)

    return

def old_gs_stuff():
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 20])
    gs.update(wspace=0, hspace=0) # set the spacing between ax
    ax = plt.subplot(gs[1])
    ax.imshow(data)
    lax = plt.subplot(gs[0])

    lclass, lcmap, lnames = convert_labels_to_cmap(row_labels, map_label_to_color)
    lax.imshow(np.matrix(lclass).T, aspect='auto', cmap=lcmap)

if __name__ == '__main__':
    data, row_colors = test_data()
    graph = sns.clustermap(data,
                           row_colors=row_colors,
                           row_cluster=False, col_cluster=False,
                           yticklabels=[''])
    _override_sns_row_colors(graph, row_colors)
    outpath = os.path.expanduser('~/tmp_sns_clustermap_hack.png')
    print "Writing to %s" %outpath
    plt.savefig(outpath, bbox_inches='tight')
    plt.close()

