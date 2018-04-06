import os
import sys
import re

import argparse
try:
    sys.path.append('/cluster/home/akahles/git/projects/2013/PanCancerTCGA/pysrc/utilities')
except IOError:
    print "Hard coded the import. I am so sorry ... "
    sys.exit()
else:
    import names
    from IPython.display import display

BASEDIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(BASEDIR)
import compute.alt_splice as preproc
import plot.alt_splice_heatmap as plotter
import config
import utils
import sf_utils

import numpy as np
import scipy.stats as stats
import pandas as pd
import h5py
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import hex2color
import seaborn as sns
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc

sys.path.append('/cluster/home/akahles/git/tools/python/viz')
import axes as axs

CACHE_DIR = os.path.expanduser('~/cache/altsplice_interest_sf')
if not os.path.exists(CACHE_DIR): os.makedirs(CACHE_DIR)

PLOT_DIR = os.path.join(config.plot_dir, 'altsplice', 'interest_sf')
if not os.path.exists(PLOT_DIR): os.makedirs(PLOT_DIR)

RESET_DF_CACHE = False

def _check_sf_interest_names(sf_names, cache_path):
    with open(cache_path, 'r') as cache_file:
        cache_names = cache_file.readlines()[0].split(',')
        cache_names = set(map(lambda x: x.lower(), cache_names))
    return cache_names == set(map(lambda x: x.lower(), sf_names))

def _cache_sf_interest_names(sf_names, cache_path):
    sf_names = sorted(map(lambda x: x.lower(), sf_names))
    with open(cache_path, 'w') as cache_file:
        cache_file.write(','.join(sf_names))
    return

def _get_interest_idx(gene_idx, ensg_ids, ensg_interest):
    gene_idx_interest = [idx for idx,name in enumerate(ensg_ids) if name.split('.')[0] in ensg_interest]
    idx_mask = np.in1d(gene_idx, gene_idx_interest)
    return np.where(idx_mask)[0]

def test_load_data(path, sf_interest):
    data = h5py.File(path, 'r')
    gene_idx = data['gene_idx'][:]
    ensg_ids = data['gene_names'][:]
    conf_idx = data['conf_idx'][:]
    interest_idx = _get_interest_idx(gene_idx, ensg_ids, sf_interest.values())
    mask_idx = np.intersect1d(conf_idx, interest_idx)
    mask = np.zeors(gene_idx.size, dtype=bool); mask[mask_idx] = True
    return data, mask

def _load_single_hdf5(path, ensg_interest):
    '''Returns the dataframe and idx of columns kept
    '''
    psi, strains, gene_idx, col_idx = sf_utils.load_data(path)
    strains = sf_utils.clean_strain(strains)
    data = h5py.File(path, 'r')
    gene_names = data['gene_names'][:]
    gene_idx_interest = [idx for idx,name in enumerate(gene_names) if name.split('.')[0] in ensg_interest]
    idx_mask = np.in1d(gene_idx, gene_idx_interest)
    psi = psi[:, idx_mask]; gene_idx = gene_idx[idx_mask]; col_idx = col_idx[idx_mask]
    columns = map(lambda x: gene_names[x.astype(int)], gene_idx)
    data.close()
    df = pd.DataFrame(psi, columns=columns, index=strains)
    return df, col_idx

def load_data_no_impute_event(etype, event_id, map_event_to_file): # , sf_interest):
    df_list = list()
    path = map_event_to_file[etype]

    data = h5py.File(path, 'r')
    gene_idx = data['gene_idx'][event_id].astype('int')
    ensg_idx = data['gene_names'][gene_idx]
    lookup = names.get_lookup_complete()
    ensg_name = names.get_ID(ensg_idx, lookup=lookup)
    psi = data['psi'][:, event_id]
    columns = [_encode_event_name(re.sub(r'-', '_', ensg_name), etype, event_id)]
    strains = sf_utils.clean_strain(data['strains'][:])
    df_list.append(pd.DataFrame(psi, index=strains, columns=columns))
    df = pd.concat(df_list, axis=1)
    return df


def load_data_no_impute(map_event_to_file, sf_interest):
    map_ensg_to_gene = {val: key for (key, val) in sf_interest.items()}
    ensg_interest = set(sf_interest.values())
    df_list = list()
    for etype, path in map_event_to_file.items():
        data = h5py.File(path, 'r')
        gene_idx = data['gene_idx'][:].astype(int)
        ensg_ids = np.vectorize(lambda x: x.split('.')[0])(data['gene_names'][:])
        ensg_idx = ensg_ids[gene_idx]
        conf_idx = data['conf_idx'][:]
        interest_idx = np.array([x for x in conf_idx if ensg_idx[x] in ensg_interest])
        mask = np.zeros(gene_idx.size, dtype=bool); mask[interest_idx] = True
        assert data['psi'].shape[1] == mask.size
        psi = data['psi'][:, mask]
        columns = [_encode_event_name(map_ensg_to_gene[ensg_idx[x]], etype, x) for x in interest_idx]
        strains = sf_utils.clean_strain(data['strains'][:])
        df_list.append(pd.DataFrame(psi, index=strains, columns=columns))
    df = pd.concat(df_list, axis=1)
    return df

def _encode_event_name(gname, etype, col_idx):
    return '%s-%s-%d' %(gname, etype, int(col_idx))

def _decode_event_name(name):
    gname, etype, col_idx = name.split('-')
    return gname, etype, int(col_idx)

def load_data(sf_interest, is_event=False):

    panc_map_event_to_file = {'exon_skip': config.alt_splice_exon_skip_path,
                             'intron_retention': config.alt_splice_intron_retention_path,
                             'alt_3prime': config.alt_splce_alt_3prime_path,
                             'alt_5prime': config.alt_splce_alt_5prime_path}

    gtex_map_event_to_file = {'exon_skip': config.alt_splice_exon_skip_gtex_path,
                             'intron_retention': config.alt_splice_intron_retention_gtex_path,
                             'alt_3prime': config.alt_splce_alt_3prime_gtex_path,
                             'alt_5prime': config.alt_splce_alt_5prime_gtex_path}


    if is_event:
        event_type = '_'.join(sf_interest.split('_')[:-1])
        event_id = int(sf_interest.split('_')[-1])
        panc_df = load_data_no_impute_event(event_type, event_id, panc_map_event_to_file)
        gtex_df = load_data_no_impute_event(event_type, event_id, gtex_map_event_to_file)
    else:
        panc_df = load_data_no_impute(panc_map_event_to_file, sf_interest)
        gtex_df = load_data_no_impute(gtex_map_event_to_file, sf_interest)
        panc_df = panc_df.loc[:, np.isnan(panc_df).mean(0) < .1]
        gtex_df = gtex_df.loc[:, np.isnan(gtex_df).mean(0) < .1]
    inter_events = panc_df.columns.intersection(gtex_df.columns)
    return panc_df.loc[:, inter_events], gtex_df.loc[:, inter_events]

def load_gtex_data(sf_interest):
    #if RESET_GTEX_CACHE: print "WARNING: Resetting gtex cache"
    sf_names = sf_interest.keys()
    ensg_interest = sf_interest.values()

    map_event_to_file = {'exon_skip': config.alt_splice_exon_skip_gtex_path,
                     'intron_retention': config.alt_splice_intron_retention_gtex_path,
                     'alt_3prime': config.alt_splce_alt_3prime_gtex_path,
                     'alt_5prime': config.alt_splce_alt_5prime_gtex_path}

    cache_path = os.path.join(CACHE_DIR, 'altsplice_interest_sf_gtex.tsv')
    cache_sf_interest_path = os.path.join(CACHE_DIR, 'altsplice_interest_sf_gtex.names')
    if False: #not RESET_GTEX_CACHE and os.path.exists(cache_path):
        assert _check_sf_interest_names(sf_names, cache_sf_interest_path)
        df = pd.read_csv(cache_path, sep='\t', index_col=0)
    else:
        # Load all data to impute nans
        # Only keep columns with genes of interest (cols describe events within gene)
        # Return intersection across all events
        df_dict = dict()
        for etype, path in map_event_to_file.items():
            df, col_idx= _load_single_hdf5(path, ensg_interest)
            df.columns = map(lambda x: etype + '.' + x, df.columns)
            df_dict[etype] = (df, col_idx)
    return df_dict

def load_panc_data(sf_interest):
    #if RESET_DF_CACHE: print "WARNING: Resetting df cache"
    sf_names = sf_interest.keys()
    ensg_interest = sf_interest.values()

    map_event_to_file = {'exon_skip': config.alt_splice_exon_skip_path,
                     'intron_retention': config.alt_splice_intron_retention_path,
                     'alt_3prime': config.alt_splce_alt_3prime_path,
                     'alt_5prime': config.alt_splce_alt_5prime_path}

    cache_path = os.path.join(CACHE_DIR, 'altsplice_interest_sf.tsv')
    cache_sf_interest_path = os.path.join(CACHE_DIR, 'altsplice_interest_sf.names')
    if False: #not RESET_DF_CACHE and os.path.exists(cache_path):
        assert _check_sf_interest_names(sf_names, cache_sf_interest_path)
        df = pd.read_csv(cache_path, sep='\t', index_col=0)
    else:
        # Load all data to impute nans
        # Only keep columns with genes of interest (cols describe events within gene)
        # Return intersection across all events
        df_dict = dict()
        for etype, path in map_event_to_file.items():
            df, col_idx = _load_single_hdf5(path, ensg_interest)
            df.columns = map(lambda x: etype + '.' + x, df.columns)
            df_dict[etype] = (df, col_idx)
    return df_dict

def filter_psi(df, min_mean, max_mean):
    means = df.mean(0)
    mask = np.logical_and(means > min_mean, means < max_mean)
    return df.loc[:, mask]

def get_sf_interest(v=False):
    gene_name_list = ['PKM', 'MAX', 'NUMB', 'MKNK2', 'BIN1',
        'RPS6KB1', 'BCL2L1', 'APC', 'PTEN',
        'KLF6', 'CD44', 'CCK2', 'GHRH',
        'CDKN2A', 'KIT']
    sf_interest = dict()
    lookup = names.get_lookup_complete()
    for gene_name in gene_name_list:
        gene_ensg = names.get_ID(gene_name, which='ensembl', lookup=lookup)
        if not gene_ensg.upper().startswith('ENSG'):
            if v: print "WARNING: %s -> %s is unexpected and ignored" %(gene_name, gene_ensg)
            continue
        sf_interest[gene_name] = gene_ensg
    if v:
        print "Genes of interest:"
        display(sf_interest)
    return sf_interest

def get_ensg_to_gene_map():
    sf_interest = get_sf_interest()
    return {val: key for key, val in sf_interest.iteritems()}

def add_col_legend(graph, col_colors):
    for label, color in sorted(col_colors.items()):
        graph.ax_col_dendrogram.bar(0, 0, color=color, label=label, linewidth=0)
    graph.ax_col_dendrogram.legend(loc='center', ncol=len(col_colors))
    return

def _get_heatmap_col_colors(psi_df):
    etype_list = ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime']
    col_cmap_lut = dict(zip(etype_list, sns.color_palette('Set2', len(etype_list)).as_hex()))
    df_etypes = psi_df.columns.map(lambda x: _decode_event_name(x)[1])
    assert np.all(map(lambda x: x in col_cmap_lut, df_etypes))
    col_colors = pd.DataFrame(df_etypes.map(lambda x: col_cmap_lut[x]), index=psi_df.columns)
    return col_colors, col_cmap_lut

def _get_heatmap_row_colors(meta_df, index):
    cnc_to_color = utils.load_color_scheme(config.color_scheme_path)
    for key, value in cnc_to_color.items():
        rgb = hex2color(value)
        cnc_to_color[key] = rgb + (1.0,)
        cnc_to_color[key + ' (Normal)'] = rgb + (.5,)
    row_colors = meta_df['cnc'].map(cnc_to_color).loc[index]
    return row_colors, cnc_to_color

def _override_sns_row_colors(graph, row_colors):

    if not isinstance(row_colors, list):
        row_colors = row_colors.tolist()
    if isinstance(row_colors[0], tuple):
        # row_colors are in rgb(a) form
        unq_colors, color_class = np.unique(row_colors, axis=0, return_inverse=True)
        unq_colors = map(lambda x: tuple(x), unq_colors)
    else:
        unq_colors, color_class = np.unique(row_colors, return_inverse=True)
        unq_colors = unq_colors.tolist()

    rcax = graph.ax_row_colors
    rcax.clear()
    cmap = colors.ListedColormap(unq_colors)
    rcax.imshow(np.matrix(color_class).T, aspect='auto', cmap=cmap)
    rcax.get_xaxis().set_visible(False)
    rcax.get_yaxis().set_visible(False)

    return


def plot_heatmap(psi_df, meta_df, outpath):
    # Sort by cancer type
    psi_df = psi_df.copy().loc[meta_df['cnc'].sort_values().index]
    psi_df = psi_df.iloc[:, psi_df.columns.map(lambda x: _decode_event_name(x)[1]).argsort()]
    col_colors, col_cmap_lut = _get_heatmap_col_colors(psi_df)
    row_colors, row_cmap_lut = _get_heatmap_row_colors(meta_df, psi_df.index)
    method = 'ward'; metric = 'cosine'
    graph = sns.clustermap(psi_df, cmap='Purples',
                           row_colors=row_colors, col_colors=col_colors,
                           row_cluster=False, col_cluster=False,
                           xticklabels=psi_df.columns.map(lambda x:_decode_event_name(x)[2]),
                           linewidths=0,
                           mask=psi_df.isnull())
    _override_sns_row_colors(graph, row_colors.values)
    graph.ax_heatmap.set_yticks([])
    graph.ax_heatmap.set_xlabel("Events")
    graph.ax_heatmap.set_ylabel("Samples")
    graph.cax.set_title("psi")
    tumor_only_row_cmap_lut = {key:val for key,val in row_cmap_lut.items() if not 'Normal' in key}
    plotter.add_legend(graph, tumor_only_row_cmap_lut)
    plotter.add_col_legend(graph, col_cmap_lut)
    print "Writing: %s" %outpath
    plt.savefig(outpath, bbox_inches='tight')
    pdf_outpath = re.sub('.png$', '.pdf', outpath)
    print "Writing: %s" %pdf_outpath
    #plt.savefig(pdf_outpath, bbox_inches='tight')
    plt.close()
    return


def run_gene_heatmaps(psi_df, mddf, sf_interest):
    plot_dir = os.path.join(PLOT_DIR, 'heatmaps')
    if not os.path.exists(plot_dir): os.makedirs(plot_dir)
    col_genes = psi_df.columns.map(lambda x: _decode_event_name(x)[0])
    for gene in sf_interest.keys():
        gene_heatmap_out = os.path.join(plot_dir, 'sf_%s_heatmap.png'%gene)
        col_mask = gene == col_genes
        if col_mask.sum() == 0: continue
        plot_heatmap(psi_df.loc[:, col_mask], mddf, gene_heatmap_out)
    return

def _get_stripplot_color_lut(unq_panc_labels, unq_gtex_labels):
    cnc_to_color = utils.load_color_scheme(config.color_scheme_path)
    for label in unq_panc_labels:
        if label in cnc_to_color: continue
        cnc = label.split()[0]
        assert cnc in cnc_to_color
        cnc_to_color[label] = cnc_to_color[cnc]
    color_lut = dict(zip(unq_gtex_labels, sns.color_palette('Set2', len(unq_gtex_labels)).as_hex()))
    color_lut.update(cnc_to_color)
    return color_lut

def run_strip_plot(panc_df, gtex_df, panc_labels, gtex_labels):
    assert panc_df.columns.equals(gtex_df.columns)
    psi_df = pd.concat((panc_df, gtex_df), axis=0)
    assert psi_df.shape[0] == panc_df.shape[0] + gtex_df.shape[0]
    assert psi_df.columns.unique().size == psi_df.shape[1]
    assert panc_df.index.equals(panc_labels.index)
    assert gtex_df.index.equals(gtex_labels.index)
    event_list = psi_df.columns.tolist()
    psi_df_aug = psi_df.copy()
    psi_df_aug['cnc'] = None
    psi_df_aug.loc[panc_df.index, ['cnc']] = panc_labels.loc[panc_df.index]
    psi_df_aug.loc[gtex_df.index, ['cnc']] = gtex_labels.loc[gtex_df.index]
    unq_panc_labels = sorted(panc_labels.unique().tolist())
    unq_gtex_labels = sorted(gtex_labels.unique().tolist())
    assert np.intersect1d(unq_panc_labels, unq_gtex_labels).size == 0
    plt.close()
    label_list = unq_panc_labels + unq_gtex_labels
    color_lut = _get_stripplot_color_lut(unq_panc_labels, unq_gtex_labels)
    for event in event_list:
        outpath = os.path.join(PLOT_DIR, 'stripplots', '%s_stripplot.png'%event)
        if not os.path.exists(os.path.dirname(outpath)): os.makedirs(os.path.dirname(outpath))
        fig, ax = plt.subplots(figsize=(20,3))
        sns.stripplot(x='cnc', y=event, data=psi_df_aug,
                      palette=color_lut, s=3,
                      order=label_list,
                      jitter=True, ax=ax)
        ax.axvline(len(unq_panc_labels) - .5, color='black', ls='--')
        ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
        ax.set_ylim(-.05, 1.05)
        ax.title.set_text('Gene: %s Event type: %s Event ID: %d'%(_decode_event_name(event)))
        ax.set_ylabel('psi')
        ax.set_xlabel('')
        axs.clean_axis(ax)
        print "Writing: %s" %outpath
        plt.savefig(outpath, bbox_inches='tight')
        pdf_outpath = re.sub('.png$', '.pdf', outpath)
        print "Writing: %s" %pdf_outpath
        plt.savefig(pdf_outpath, bbox_inches='tight')
        plt.close()
    return

def filter_both_psi(panc_df, gtex_df, min_mean=.05, max_mean=.95):
    panc_mean = stats.nanmean(panc_df.values, 0)
    gtex_mean = stats.nanmean(gtex_df.values, 0)
    mask = (panc_mean > min_mean) * (panc_mean < max_mean) * (gtex_mean > min_mean) * (gtex_mean < max_mean)
    return panc_df.loc[:, mask], gtex_df.loc[:, mask]

def load_gtex_metadata(path, index=None):
    df = pd.read_csv(path, sep='\t')[['Run_s', 'histological_type_s']].set_index('Run_s')
    if index is not None:
        df = df.loc[index]
        assert df.index.size == index.size
    return df

def add_normal_label(md):
    normal_mask = np.invert(md['is_tumor'])
    md.loc[normal_mask, 'cnc'] = md.loc[normal_mask]['cnc'].map(lambda x: x+' (Normal)')
    return md

DEBUG = False
if __name__ == '__main__':

    if sys.argv[1] == '--event':

        ### load actual data
        panc_df, gtex_df = load_data(sys.argv[2], is_event=True)

        ### load metadata
        panc_md = utils.load_metadata_df(config.metadata_path, panc_df.index).dropna()
        gtex_md = load_gtex_metadata(config.gtex_metadata_path, gtex_df.index).dropna()
        panc_md_wnorm = add_normal_label(panc_md)
        assert panc_df.index.equals(panc_md.index)
        assert gtex_df.index.equals(gtex_md.index)

        if len(sys.argv) > 3:
            PLOT_DIR = os.path.join(config.plot_dir, 'altsplice', 'outlier_events_%s' % sys.argv[3])
        else:    
            PLOT_DIR = os.path.join(config.plot_dir, 'altsplice', 'outlier_events')
        if not os.path.exists(PLOT_DIR): os.makedirs(PLOT_DIR)

        run_strip_plot(panc_df, gtex_df, panc_md_wnorm['cnc'], gtex_md['histological_type_s'])
    else:
        sf_interest = get_sf_interest()

        ### load actual data
        panc_df, gtex_df = load_data(sf_interest)

        ### load metadata
        panc_md = utils.load_metadata_df(config.metadata_path, panc_df.index).dropna()
        gtex_md = load_gtex_metadata(config.gtex_metadata_path, gtex_df.index).dropna()
        panc_md_wnorm = add_normal_label(panc_md)
        assert panc_df.index.equals(panc_md.index)
        assert gtex_df.index.equals(gtex_md.index)

        run_gene_heatmaps(panc_df, panc_md, sf_interest)
