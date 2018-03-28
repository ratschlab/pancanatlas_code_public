import os
import sys
import pandas as pd
import numpy as np
import h5py
import re

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import sqtl
import alt_splice_heatmap_hl as hmhelp
import alt_splice_embeddings as ebhelp

BASEDIR = os.path.dirname(os.path.dirname(__file__))
sys.path.append(BASEDIR)
import compute.alt_splice as preproc
import config
import utils

import scipy.spatial as sp
import scipy.cluster.hierarchy as hc

def get_ensg_of_interest(genes_with_targets, trans_df):
    return np.asarray(trans_df[trans_df.loc[:, 'gene-id(mutation)'].isin(genes_with_targets)].index.unique())

def get_event_pos_of_interest(genes_with_targets, trans_df):
    map_gwt_to_event_pos = dict()
    event_pos_cols = [col for col in trans_df.columns if 'event-pos' in col]
    for gwt in genes_with_targets:
        mask = trans_df.loc[:, 'gene-id(mutation)'] == gwt
        gwt_event_pos_df = trans_df[mask].loc[:, event_pos_cols]
        gwt_event_pos = np.unique(['_'.join(x) for x in gwt_event_pos_df.astype(str).values])
        map_gwt_to_event_pos[gwt] = gwt_event_pos
    return map_gwt_to_event_pos

def load_focused_events(gtdf, ensg_list):

    df_cache = os.path.expanduser('~/cache/alt_splice_heatmap/sqtl/heatmap_focused.tsv')
    ensg_map_cache = os.path.expanduser('~/cache/alt_splice_heatmap/sqtl/heatmap_focused_ensg_map.npy')

    if False and os.path.exists(df_cache) and os.path.exists(ensg_map_cache):
        df = pd.read_csv(df_cache, sep='\t')
        ensg_map = np.load(ensg_map_cache).item()
    else:
        map_etype_to_file = {'exon_skip': config.alt_splice_exon_skip_path,
#                             'intron_retention': config.alt_splice_intron_retention_path,
                             'alt_3prime': config.alt_splce_alt_3prime_path,
                             'alt_5prime': config.alt_splce_alt_5prime_path}

        df_list = list()
        ensg_map = dict()
        for event, path in map_etype_to_file.items():
            print "Loading %s from %s" %(event, path)
            psi, strains, gene_idx = preproc.load_data(path)
            gene_idx = gene_idx.astype(int)
            gene_names = h5py.File(path, 'r')['gene_names'][:]
            gi_mask = np.in1d(map(lambda x: x.split('.')[0], gene_names[gene_idx]), ensg_list)
            gene_idx = gene_idx[gi_mask]
            psi = psi[:, gi_mask]
            columns = [event + '-' + str(int(x)) for x in gene_idx]
            ensg_map.update(dict(zip(columns, gene_names[gene_idx])))
            strains = preproc.clean_strain(strains)
            df = pd.DataFrame(psi, index=strains, columns=columns)
            df_list.append(df)
        df = pd.concat(df_list, axis=1, join='inner')
        if not os.path.exists(os.path.dirname(df_cache)): os.makedirs(os.path.dirname(df_cache))
        df.to_csv(df_cache, sep='\t')
        np.save(ensg_map_cache, ensg_map)

    return df, ensg_map

def do_heatmap(gtdf, psi_df, map_gwt_to_event_pos, column_events):
    column_events = np.array(column_events)
    metadata_df = utils.load_metadata_df(config.metadata_path, psi_df.index)
    gtdf = sqtl._match_gtid_to_data_index(gtdf, psi_df.index)

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
    column_events = column_events[low_var_mask]

    # Get color scheme
    col_keys = np.unique(psi_df.columns.map(lambda x: x.split('-')[0]))
    cnc_to_color = utils.load_color_scheme(config.color_scheme_path)
    col_cmap = np.array(sns.color_palette('Set2', len(col_keys)))
    col_cmap_lut = dict(zip(col_keys, col_cmap))

    nan_class = 3
    if np.any(gtdf.isnull()):
        assert not nan_class in gtdf.values
        gtdf.values[np.isnan(gtdf.values)] = nan_class
    gtdf = gtdf.apply(lambda x: pd.to_numeric(x, downcast='integer'))
    all_col_colors = hmhelp.ash.map_col_colors(psi_df.columns, col_cmap_lut)

    TEST_MODE = False
    if TEST_MODE:
        print "WARNING: TEST MODE!"

    ensg_list = ['ENSG00000115524', 'ENSG00000138413', 'ENSG00000054282']
    half_height = True
    for mask_to_events in [False]:
#        for norm_method in ['none', 'meanzero', 'standardize']:
         for norm_method in ['standardize']:
            #genes_adhoc = [u'ENSG00000115524.11_2_198266834', u'ENSG00000115524.11_2_198267483', u'ENSG00000138413.9_2_209113112', u'ENSG00000054282.11_1_243493891']
            #for gene in genes_adhoc:
            for gene in gtdf.columns:
                ensg, chrm, pos = gene.split('_')
                ensg = ensg.split('.')[0]
                if not (ensg in ensg_list): continue
                if mask_to_events:
                    gene_event_pos = map_gwt_to_event_pos[gene.split('_')[0].replace('-','/')]
                    gene_psi_df = psi_df.iloc[:, np.in1d(column_events, gene_event_pos)]
                else:
                    gene_psi_df = psi_df.copy()

                # sort by mutation status
                mut_status = gtdf[gene][gtdf[gene] != nan_class].sort_values()
                unq_mut_status = mut_status.unique()
                gene_psi_df = gene_psi_df.loc[mut_status.index]
                row_cscheme = sns.color_palette('Set2', unq_mut_status.size + len(col_keys))[len(col_keys):]
                row_cmap = dict(zip(unq_mut_status, row_cscheme))
                all_row_colors = mut_status.map(lambda x:row_cmap[x])
                row_labels = ['Alternate', 'Heterozygous', 'Reference']
                map_label_to_color = dict([(row_labels[x], row_cmap[x]) for x in unq_mut_status])
                #for viz in ['full', 'sampled', 'averaged']:
                for viz in ['sampled']:

                    psi_df_desc = 'focused_%s' %viz
                    if mask_to_events:
                        psi_df_desc = psi_df_desc + '_' + 'masked'
                    outpath = os.path.join(_HEATMAP_BASE, ensg, gene.replace('.','-') + '_' + psi_df_desc + '.png')

                    if viz == 'full':
                        graph_df = gene_psi_df.loc[mut_status.index]

                    counts = mut_status.value_counts()
                    if viz == 'sampled':
                        min_count = counts.iloc[counts.index != nan_class].min()
                        if min_count < 50: min_count = 50
                        index = list()
                        npop = dict()
                        for val, cnt in counts.sort_index().iteritems():
                            mask = mut_status == val
                            npop[val] = mask.sum()
                            nsamples = min(min_count, mask.sum())
                            index.extend(np.random.choice(mut_status.index[mask], nsamples, replace=False))
                        psi_df_desc = psi_df_desc + '_nmutated_%d_nref_%d' %(npop[1], npop[2])
                        graph_df = gene_psi_df.loc[index]

                    if viz == 'averaged':
                        avgs = list()
                        index = list()
                        for val, cnt in counts.sort_index().iteritems():
                            mask = mut_status == val
                            if mask.sum() == 0: continue
                            avgs.append(gene_psi_df.loc[mask].mean(0))
                            index.append(gene_psi_df.loc[mask].index[0])
                        graph_df = pd.concat(avgs, axis=1).T
                        graph_df.index = index

                    cax_title = 'psi'
                    if norm_method == 'standardize':
                        norm = (graph_df.values - graph_df.values.mean(0)) / graph_df.values.std(0)
                        limit = 1
                        norm[np.isnan(norm)] = 0
                        graph_df.iloc[:] = norm
                        psi_df_desc = psi_df_desc + '_standardized'
                        cax_title = 'Z-score (PSI)'

                    elif norm_method == 'meanzero':
                        norm = graph_df.values - graph_df.values.mean(0)
                        norm[np.isnan(norm)] = 0
                        graph_df.iloc[:] = norm
                        psi_df_desc = psi_df_desc + '_meanzero'
                        cax_title = 'meanzero psi'

                    col_linkage = hc.linkage(sp.distance.pdist(graph_df.values.T), method='ward', metric='cosine')
                    row_colors = all_row_colors.loc[graph_df.index]
                    col_colors =  all_col_colors.loc[graph_df.columns]
                    col_colors.name = None
                    row_colors.name = None

                    # And finally plot the data
                    print "Plotting data ... "
                    print "graph_df.shape = " + str(graph_df.shape)
                    figsize = None
                    if half_height:
                        figsize = (10, 5)
                        psi_df_desc += '_halfheight'
                    graph = sns.clustermap(graph_df,
                                       row_colors=row_colors.values.tolist(), col_colors=col_colors.values.tolist(),
                                       row_cluster=False, col_linkage=col_linkage,
                                       cbar_kws={"ticks":[-1, 0, 1]},
                                       vmin=-limit, vmax=limit,
                                       cmap='BrBG',
                                       figsize=figsize)

                    graph.ax_heatmap.xaxis.set_ticklabels([])
                    graph.ax_heatmap.yaxis.set_ticklabels([])
                    graph.ax_heatmap.xaxis.set_ticks([])
                    graph.ax_heatmap.yaxis.set_ticks([])

                    # Remove column dendrogram
                    graph.ax_col_dendrogram.set_xlim([0,0])

                    mut_event = "Gene: %s, Chr %s Loc %s" %(ensg, chrm, pos)
                    #graph.ax_col_dendrogram.set_title("%s AltSplice %s Clustering" %(mut_event, psi_df_desc.replace('_', ' ').title()))
                    graph.ax_heatmap.set_xlabel("Events")
                    graph.ax_heatmap.set_ylabel("Samples")
                    graph.cax.set_title(cax_title)
                    if TEST_MODE or ensg == ensg_list[-1]:
                        # add legends
                        hmhelp.ash.add_legend(graph, map_label_to_color)
                        hmhelp.ash.add_col_legend(graph, col_cmap_lut)
                    graph.ax_col_colors.set_yticks([])
                    graph.ax_row_colors.set_xticks([])
                    graph.ax_row_colors.set_ylabel('chr%s %s'%(chrm, pos))

                    #Removing colorbar!
                    if TEST_MODE or ensg == ensg_list[-1]:
                        graph.cax.set_yticklabels(['< -1', '0', '> 1'])
                    else:
                        graph.cax.set_visible(False)

                    outpath = os.path.join(_HEATMAP_BASE, ensg, gene.replace('.','-') + '_' + psi_df_desc +  '_v3.png')
                    if TEST_MODE:
                        print "WARNING: TEST MODE ACTIVE"
                        outpath = os.path.join(os.path.expanduser('~/tmp'), os.path.basename(outpath))
                    if not os.path.exists(os.path.dirname(outpath)): os.makedirs(os.path.dirname(outpath))

                    print "Saving heatmap to: %s" % outpath
                    print "CONTINUE"
                    continue
                    plt.savefig(outpath, bbox_inches='tight', dpi=500)
                    print "Saving heatmap to: %s" % re.sub(r'png$', 'pdf', outpath)
                    plt.savefig(re.sub(r'png$', 'pdf', outpath), format='pdf', bbox_inches='tight')
                    plt.close()

_HEATMAP_BASE = os.path.join(config.plot_dir, 'altsplice', 'concatenated', 'heatmaps', 'sqtl')

def load_psi_df(map_gwt_to_event_pos):
    event_pos = np.unique([item for sublist in map_gwt_to_event_pos.values() for item in sublist])
    # load data, impute values, filter to event pos
    kjong_psi_h5_path = '/cluster/work/grlab/projects/TCGA/PanCancer/hdf5_10Kset/phenotypes.hdf5r10_s2000_V100.hdf5'
    psih5 = h5py.File(kjong_psi_h5_path, 'r')

    column_events = list()
    df_cache = os.path.expanduser('~/cache/alt_splice_heatmap/sqtl/heatmap_focused.tsv')
    col_events_cache = os.path.expanduser('~/cache/alt_splice_heatmap/sqtl/col_events.npy')
    if True or not os.path.exists(df_cache):
        df_list = list()
        for key in  ['exon_skip', 'alt_3prime', 'alt_5prime']:
            key_event_pos = ['_'.join(x) for x in psih5['splicing'][key]['event_pos'][:].T.astype(str)]
            mask = np.in1d(key_event_pos, event_pos)
            if mask.sum() == 0: continue
            assert mask.size == psih5['splicing'][key]['psi'].shape[0]
            psi = psih5['splicing'][key]['psi'][mask, :]
            strains = psih5['splicing'][key]['strains'][:]
            gene_idx = psih5['splicing'][key]['gene_idx'][mask]
            used_key_event_pos = np.array(key_event_pos)[mask]

            psi, strains, gene_idx, col_idx = preproc.preproc(psi.T, strains, gene_idx, ret_col_idx=True)
            strains = preproc.clean_strain(strains)
            assert np.all(np.in1d(used_key_event_pos[col_idx], event_pos))
            column_events.extend(used_key_event_pos[col_idx])
            columns = [key + '-' + str(int(x)) for x in gene_idx]
            df_list.append(pd.DataFrame(psi, index=strains, columns=columns))
            assert psi.shape[1] == col_idx.size
        psi_df = pd.concat(df_list, axis=1, join='inner')
        assert psi_df.shape[1] == len(column_events)
        psi_df.to_csv(df_cache, sep='\t')
        np.save(col_events_cache, column_events)
    else:
        print "Loading psi_df from cache: %s" %df_cache
        psi_df = pd.read_csv(df_cache, sep='\t', index_col=0)
        column_events = np.load(col_events_cache)
    return psi_df, column_events

def load_trans_info(trans_path=None):
    if trans_path is None:
        trans_path = sqtl._TRANS_DATA_PATH
    trans_df = pd.read_csv(trans_path, sep='\t')
    if trans_df.columns[0].startswith('#'):
        cols = trans_df.columns.tolist()
        cols[0] = cols[0].replace('# ', '')
        trans_df.columns = cols
    trans_df = trans_df.set_index('gene-id(event)')

    # find all genes with more than k=50 downstream targets
    genes_with_targets = sqtl.find_genes_with_targets(trans_df, k=50)
    event_pos = get_event_pos_of_interest(genes_with_targets, trans_df)
    return trans_df, genes_with_targets, event_pos

def load_mutation_info(trans_df, genes_with_targets):
    # get all coords
    # for all coords get mutation information for all samples
    coords = sqtl.get_coords(trans_df, genes_with_targets)
    gt_file = h5py.File(sqtl._GENOTYPE_H5_PATH, 'r')
    gtdf = sqtl.get_mutation_info(gt_file, coords)
    return gtdf

if __name__ == '__main__':
    trans_df, genes_with_targets, map_gwt_to_event_pos = load_trans_info()
    gtdf = load_mutation_info(trans_df, genes_with_targets)
    psi_df, column_events = load_psi_df(map_gwt_to_event_pos)
    do_heatmap(gtdf, psi_df, map_gwt_to_event_pos, column_events)

