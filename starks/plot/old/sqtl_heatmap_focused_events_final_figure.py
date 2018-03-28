import os
import sys
import pandas as pd
import numpy as np
import re
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

import sqtl_heatmap_focused_events as helper
config = helper.config
utils = helper.utils

sys.path.append('/cluster/home/starks/git/tools-python/viz')
import heatmap

MUT_STATUS_NAN_CLASS = 3

def _format_ax(ax):
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    return

def _remove_spines(ax):
    for loc in ['top', 'right', 'bottom', 'left']:
        ax.spines[loc].set_visible(False)
    return

def load_and_preproc():
    trans_df, genes_with_targets, map_gwt_to_event_pos = helper.load_trans_info()
    gtdf = helper.load_mutation_info(trans_df, genes_with_targets)
    psi_df, column_events = helper.load_psi_df(map_gwt_to_event_pos)

    column_events = np.array(column_events)
    metadata_df = helper.utils.load_metadata_df(helper.config.metadata_path, psi_df.index)
    gtdf = helper.sqtl._match_gtid_to_data_index(gtdf, psi_df.index)

    # Drop missing samples
    index_with_gt_data = gtdf.index[~np.all(gtdf.isnull(), axis=1)]
    index_with_all_data = index_with_gt_data.intersection(metadata_df.index)
    psi_df = psi_df.loc[index_with_all_data]
    gtdf = gtdf.loc[index_with_all_data]
    metadata_df = metadata_df.loc[index_with_all_data]

    # Order by cancer type, event_type
    psi_df = psi_df.loc[metadata_df['cnc'].sort_values().index]
    psi_df = psi_df.iloc[:, psi_df.columns.argsort()]
    column_events = column_events[psi_df.columns.argsort()]

    # Remove low var events
    low_var_mask = psi_df.std() > .01
    psi_df = psi_df.loc[:, low_var_mask]

    if np.any(gtdf.isnull()):
        assert not MUT_STATUS_NAN_CLASS in gtdf.values
        gtdf.values[np.isnan(gtdf.values)] = MUT_STATUS_NAN_CLASS
    gtdf = gtdf.apply(lambda x: pd.to_numeric(x, downcast='integer'))

    return psi_df, gtdf

def _get_col_keys(psi_df):
    '''extract alt_3prime, alt_5prime, exon_skip, etc.'''
    col_keys = psi_df.columns.map(lambda x: x.split('-')[0])
    return col_keys

def load_cnc_colors(psi_df):
    # Get color scheme
    col_keys = _get_col_keys(psi_df)
    col_cmap = np.array(sns.color_palette('Set2', len(col_keys)))
    col_cmap_lut = dict(zip(col_keys, col_cmap))

def single_heatmap(psi_df, mut_status, axes):
    ax_hmap, ax_row_color, ax_col_color = axes
    # sort by mutation status
    psi_df = psi_df.loc[mut_status.index]

    # sample based on mut_status
    viz = 'sampled'
    psi_df_desc = 'focused_%s' %viz
    counts = mut_status.value_counts()
    min_count = counts.iloc[counts.index != MUT_STATUS_NAN_CLASS].min()
    if min_count < 50: min_count = 50
    index = list()
    nsample_dict = dict()
    for val, cnt in counts.sort_index().iteritems():
        mask = mut_status == val
        nsamples = min(min_count, mask.sum())
        nsample_dict[val] = nsamples
        index.extend(np.random.choice(mut_status.index[mask], nsamples, replace=False))
    psi_df = psi_df.loc[index]

    # Standardize
    cax_title = 'psi'
    norm = (psi_df.values - psi_df.values.mean(0)) / psi_df.values.std(0)
    norm[np.isnan(norm)] = 0
    psi_df.iloc[:] = norm
    psi_df_desc = psi_df_desc + '_standardized'
    cax_title = 'Z-score (PSI)'

    # Cluster columns
    col_linkage = hc.linkage(sp.distance.pdist(psi_df.values.T), method='ward', metric='cosine')
    dendro = hc.dendrogram(col_linkage, p=100000, no_plot=True, truncate_mode='mtica')
    lvs = dendro['leaves']
    psi_df = psi_df.iloc[:, lvs]

    # color cols from event type
    event_type = _get_col_keys(psi_df)
    et_vals, et_classes = np.unique(event_type, return_inverse=True)
    et_cmap = ListedColormap(sns.color_palette('Set2', et_vals.size))
    et_lut = dict([(et, et_cmap(i)) for i,et in enumerate(et_vals)])

    # color rows from mut_status
    mut_labels = ['Alternate', 'Heterozygous', 'Reference', 'UNK']
    mut_vals, mut_classes = np.unique(mut_status.loc[psi_df.index], return_inverse=True)
    mut_cmap = ListedColormap(sns.color_palette('Set2', mut_vals.size + et_vals.size)[et_vals.size:])
    mut_lut = dict([(mut_labels[indx], mut_cmap(i)) for i,indx in enumerate(mut_vals)])

    print "Plotting data ... "
    print "psi_df.shape = " + str(psi_df.shape)

    hmap = ax_hmap.imshow(psi_df.values, aspect='auto', cmap='BrBG', vmin=-1.0, vmax=1.0)
    ax_col_color.imshow(np.matrix(et_classes), aspect='auto', cmap=et_cmap)
    ax_row_color.imshow(np.matrix(mut_classes).T, aspect='auto', cmap=mut_cmap)
    return hmap, et_lut, mut_lut

_OUTDIR = os.path.join(config.plot_dir,
                'altsplice', 'concatenated', 'heatmaps', 'sqtl', 'final_figs')

def main_fig_proc(psi_df, gtdf):
    fig = plt.figure(figsize=(10, 16))
    outer_grid = gridspec.GridSpec(5, 2,
                                   width_ratios=[1, 30],
                                   height_ratios=[7, 7, 7, .5, .5],
                                   hspace=0.25, wspace=0.1)

    ensg_pos_list = ['ENSG00000115524.11_2_198267483',
                     'ENSG00000138413.9_2_209113112',
                     'ENSG00000054282.11_1_243493891']

    for i, ensg_pos in enumerate(ensg_pos_list):
        ensg, chrm, pos = ensg_pos.split('_')
        inner_grid = gridspec.GridSpecFromSubplotSpec(2, 2,
                width_ratios=[1,20],
                height_ratios=[1,20],
                subplot_spec=outer_grid[i, 1:],
                wspace=0.025, hspace=0.05)

        ax_col_color = plt.subplot(inner_grid[0,1]); _format_ax(ax_col_color); _remove_spines(ax_col_color)
        ax_hmap = plt.subplot(inner_grid[1,1]); _format_ax(ax_hmap)
        ax_row_color = plt.subplot(inner_grid[1,0]); _format_ax(ax_row_color); _remove_spines(ax_row_color)
        axes = (ax_hmap, ax_row_color, ax_col_color)

        mut_status = gtdf[ensg_pos][gtdf[ensg_pos] != MUT_STATUS_NAN_CLASS].sort_values()

        hmap, col_lut, row_lut = single_heatmap(psi_df.copy(), mut_status, axes)
        ax_row_color.set_ylabel('chr%s %s'%(chrm, pos), fontsize=10)
        ax_hmap.set_xlabel('Events', fontsize=15)
        ax_hmap.yaxis.set_label_position("right")
        ax_hmap.set_ylabel('Samples', fontsize=15)

    ax_cbar = plt.subplot(outer_grid[-2:, 0]); _format_ax(ax_cbar); _remove_spines(ax_cbar)
    ax_xleg = plt.subplot(outer_grid[-2, 1]);  _format_ax(ax_xleg);
    ax_yleg = plt.subplot(outer_grid[-1, 1]);  _format_ax(ax_yleg);

    fs = 10
    fig.colorbar(hmap, cax=ax_cbar, ticks=[-1,0,1])
    ax_cbar.set_yticklabels(['< -1', '0', '> 1'], fontsize=fs)
    ax_cbar.set_title('Z-score (PSI)', fontsize=fs)

    fs = 15
    col_lgd_patches = [mpatches.Patch(color=clr, label=lbl) for lbl, clr in col_lut.items()]
    ax_xleg.axis('off')
    ax_xleg.legend(handles=col_lgd_patches, ncol=len(col_lut), loc='upper center', fontsize=fs)

    row_lgd_patches = [mpatches.Patch(color=clr, label=lbl) for lbl, clr in row_lut.items()]
    ax_yleg.axis('off')
    ax_yleg.legend(handles=row_lgd_patches, ncol=len(row_lut), loc='upper center', fontsize=fs)

    outpath = os.path.join(_OUTDIR, 'fig4.png')
    if not os.path.exists(os.path.dirname(outpath)): os.makedirs(os.path.dirname(outpath))

    print "Saving heatmap to: %s" % outpath
    plt.savefig(outpath, bbox_inches='tight', dpi=200)
#    print "Saving heatmap to: %s" % re.sub(r'png$', 'pdf', outpath)
#    plt.savefig(re.sub(r'png$', 'pdf', outpath), format='pdf', bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    psi_df, gtdf = load_and_preproc()
#    main_fig_proc(psi_df, gtdf)

