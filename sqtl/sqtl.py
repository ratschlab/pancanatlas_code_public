import os
import sys
import pandas as pd
import numpy as np
import h5py

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import alt_splice_heatmap_hl as hmhelp
import alt_splice_embeddings as ebhelp

BASEDIR = os.path.dirname(os.path.dirname(__file__))
sys.path.append(BASEDIR)
import compute.alt_splice as preproc
import config
import utils

import scipy.spatial as sp
import scipy.cluster.hierarchy as hc

_HEATMAP_BASE = os.path.join(config.plot_dir, 'altsplice', 'concatenated', 'heatmaps', 'sqtl')
_EMBED_BASE = os.path.join(config.plot_dir, 'altsplice', 'concatenated', 'embeds', 'sqtl')

#_TRANS_DATA_PATH = '/cluster/work/grlab/projects/TCGA/PanCancer/QTL_Analysis_results/analysis_10kset/transdataresults.tsv'
_TRANS_DATA_PATH = '/cluster/work/grlab/projects/TCGA/PanCancer/QTL_Analysis_results/analysis_10kset/2dmanhattan_information_v3/transdataresults_test_v3.tsv'
_GENOTYPE_H5_PATH  = '/cluster/work/grlab/projects/TCGA/PanCancer/hdf5_10Kset/mergedFiles_clean.hdf5'

def find_genes_with_targets(trans_df, k=50):
    targets = trans_df.loc[:, 'gene-id(mutation)']
    genes, counts = np.unique(targets.values, return_counts=True)
    return genes[counts >= k]

def get_coords(trans_df, target_genes):
    target_mask = trans_df.loc[:, 'gene-id(mutation)'].isin(target_genes)
    coords = trans_df.loc[target_mask, ['snp_chrm', 'snp_pos', 'gene-id(mutation)']]
    coords = coords.drop_duplicates()
    coords.index = np.arange(coords.shape[0])
    coords.columns = ['chr', 'pos', 'gene']
    return coords

def _name_coords(coords):
    names = list()
    for indx, row in coords.iterrows():
        names.append('_'.join([row['gene'].replace('/', '-'), str(row['chr']), str(row['pos'])]))
    return names

def get_mutation_info(gt_file, coords):
    '''DataFram of gt data for positions located in coords (index=samples, cols=locs).
    '''
    gt_pos = gt_file['pos'][:]
    use_gt_pos = list()
    for (gchr, gpos) in coords[['chr', 'pos']].values:
        mask = (gt_pos[:, 0] == gchr) * (gt_pos[:, 1] == gpos)
        assert mask.sum() == 1
        use_gt_pos.append(np.where(mask)[0][0])
    use_gt_pos = np.array(use_gt_pos)
    coords = coords.iloc[np.argsort(use_gt_pos)]
    names = _name_coords(coords)
    mask = np.zeros(gt_file['pos'].shape[0], dtype=bool)
    mask[use_gt_pos] = True
    data = gt_file['gt'][mask, :].T
    gtids = gt_file['gtid'][:]
    return pd.DataFrame(data, index=gtids, columns=names)

def load_gt_data(k=50):
    trans_df = pd.read_csv(_TRANS_DATA_PATH, sep='\t')
    if trans_df.columns[0].startswith('#'):
        cols = trans_df.columns.tolist()
        cols[0] = cols[0].replace('# ', '')
        trans_df.columns = cols
    trans_df = trans_df.set_index('gene-id(event)')

    # find all genes with more than k=50 downstream targets
    genes_with_targets = find_genes_with_targets(trans_df, k=k)

    # get all coords
    coords = get_coords(trans_df, genes_with_targets)

    # for all coords get mutation information for all samples
    gt_file = h5py.File(_GENOTYPE_H5_PATH, 'r')
    gt_df = get_mutation_info(gt_file, coords)

    return gt_df


def _match_gtid_to_data_index(gtdf, index):
    gtdf = gtdf.copy()
    gtdf.index = gtdf.index.map(lambda x: '-'.join(x.split('-')[:3]))
    gtdf = gtdf.loc[index.map(lambda x: '-'.join(x.split('-')[:3]))]
    gtdf.index = index
    return gtdf

def map_col_colors(event_names, lut):
    col_colors = [lut[ev.split('-')[0]] for ev in event_names]
    col_colors = pd.Series(col_colors, index=event_names)
    return col_colors

def get_col_linkage(combined_df, method='ward', metric='cosine'):
    CACHE_DIR = os.path.expanduser('~/cache/alt_splice_heatmap/sqtl')
    if not os.path.exists(CACHE_DIR): os.makedirs(CACHE_DIR)
    col_linkage_cache_path = os.path.join(CACHE_DIR, 'col_linkage_%s_%s.npy' %(method, metric))
    idx_linkage_cache_path = os.path.join(CACHE_DIR, 'idx.npy')
    col_name_cache_path = os.path.join(CACHE_DIR, 'col_names.npy')
    if os.path.exists(col_linkage_cache_path):
        print "Loading linkage from %s" %col_linkage_cache_path
        col_linkage = np.load(col_linkage_cache_path)
        assert np.array_equal(np.load(idx_linkage_cache_path), combined_df.index)
        assert np.array_equal(np.load(col_name_cache_path), combined_df.columns)
    else:
        print "Calculating linkage"
        col_linkage = hc.linkage(sp.distance.pdist(combined_df.values.T), method=method, metric=metric)
        np.save(col_linkage_cache_path, col_linkage)
        np.save(idx_linkage_cache_path, combined_df.index)
        np.save(col_name_cache_path, combined_df.columns)
    return col_linkage

def do_heatmap_full(gtdf):

    # Load the data
    combined_df_desc = '%s_%d_high_var_events_concat' %(etype_desc, MAX_EVENTS)
    combined_df = hmhelp.load_full_dataset(map_etype_to_file, MAX_EVENTS, combined_df_desc, reset=False)

    metadata_df = utils.load_metadata_df(config.metadata_path, combined_df.index)
    gtdf = _match_gtid_to_data_index(gtdf, combined_df.index)

    # Drop missing samples
    index_with_gt_data = gtdf.index[~np.all(gtdf.isnull(), axis=1)]
    index_with_all_data = index_with_gt_data.intersection(metadata_df.index)
    combined_df = combined_df.loc[index_with_all_data]
    gtdf = gtdf.loc[index_with_all_data]
    metadata_df = metadata_df.loc[index_with_all_data]

    # Order by cancer type, event_type
    combined_df = combined_df.loc[metadata_df['cnc'].sort_values().index]
    combined_df = combined_df.iloc[:, combined_df.columns.argsort()]

    # Get color scheme
    cnc_to_color = utils.load_color_scheme(config.color_scheme_path)
    col_cmap = np.array(sns.color_palette('Set2', len(map_etype_to_file.keys())))
    col_cmap_lut = dict(zip(map_etype_to_file.keys(), col_cmap))


    gtdf = gtdf.copy()
    nan_class = 3
    if np.any(gtdf.isnull()):
        assert not nan_class in gtdf.values
        gtdf.values[np.isnan(gtdf.values)] = nan_class
    gtdf = gtdf.apply(lambda x: pd.to_numeric(x, downcast='integer'))
    row_cmap = sns.color_palette('Set2', nan_class+1)
    all_col_colors = hmhelp.ash.map_col_colors(combined_df.columns, col_cmap_lut)
    full_col_linkage = get_col_linkage(combined_df)


    for gene in gtdf.columns:
        pass
    for gene in ['ENSG00000138413.9_2_209113112']:
        # TODO: normalize -- subtract mean and divide by stdv
        # TODO: change color
        # sort by mutation status
        mut_status = gtdf[gene].sort_values()
        all_row_colors = mut_status.map(lambda x:row_cmap[x])

        for viz in ['full', 'sampled', 'averaged']:
            if viz == 'full':
                continue;
                graph_df = combined_df.loc[mut_status.index]

            counts = mut_status.value_counts()
            if viz == 'sampled':
                min_count = counts.iloc[counts.index != nan_class].min()
                if min_count < 50: min_count = 50
                index = list()
                for val, cnt in counts.sort_index().iteritems():
                    mask = mut_status == val
                    nsamples = min(min_count, mask.sum())
                    index.extend(np.random.choice(mut_status.index[mask], nsamples, replace=False))
                graph_df = combined_df.loc[index]

            if viz == 'averaged':
                avgs = list()
                index = list()
                for val, cnt in counts.sort_index().iteritems():
                    mask = mut_status == val
                    avgs.append(combined_df.loc[mask].mean(0))
                    index.append(combined_df.loc[mask].index[0])
                graph_df = pd.concat(avgs, axis=1).T
                graph_df.index = index

            graph_df.iloc[:] = (graph_df.values - graph_df.values.mean(0) )
            row_colors = all_row_colors.loc[graph_df.index]
            col_colors =  all_col_colors.loc[graph_df.columns]

            # And finally plot the data
            sys.setrecursionlimit(100000)
            print "Plotting data ... "
            graph = sns.clustermap(graph_df,
                               row_colors=row_colors, col_colors=col_colors,
                               row_cluster=False, col_linkage=full_col_linkage,
                               cmap='BrBG')

            graph.ax_heatmap.xaxis.set_ticklabels([])
            graph.ax_heatmap.yaxis.set_ticklabels([])
            graph.ax_heatmap.xaxis.set_ticks([])
            graph.ax_heatmap.yaxis.set_ticks([])


            mut_event = "Gene: %s, Chr %s Loc %s" %tuple(gene.split('_'))
            graph.ax_col_dendrogram.set_title("%s AltSplice %s Clustering" %(mut_event, combined_df_desc.replace('_', ' ').title()))
            graph.ax_heatmap.set_xlabel("Events")
            graph.ax_heatmap.set_ylabel("Samples")
            graph.cax.set_title("psi")
            row_labels = map(str, range(nan_class)) + ['Missing']
            hmhelp.ash.add_legend(graph, dict(zip(row_labels, row_cmap)))
            hmhelp.ash.add_col_legend(graph, col_cmap_lut)

            ensg = gene.split('_')[0].split('.')[0]
            outpath = os.path.join(_HEATMAP_BASE, ensg, gene.replace('.','-') + '_' + combined_df_desc + '_v%d.png'%_VERSION)
            if viz != 'full':
                outpath = os.path.join(os.path.dirname(outpath), viz + '_' + os.path.basename(outpath))
            if not os.path.exists(os.path.dirname(outpath)): os.makedirs(os.path.dirname(outpath))

            print "Saving heatmap to: %s" %outpath
            plt.savefig(outpath, bbox_inches='tight', dpi=300)
            plt.close()
            if DEBUG and viz == 'sampled': return

def _copy_and_update(copy_dict, update_dict):
    ret_dict = copy_dict.copy()
    ret_dict.update(update_dict)
    return ret_dict

def do_embedding(gtdf):
    gtdf = gtdf.copy()
    DO_TSNE = True

    DO_INDIVIDUAL_EVENTS = True
    DO_COMBINED_EVENTS = True

    PLOT_DIR = os.path.join(config.plot_dir, 'altsplice')
    EMBED_DIR = os.path.join(config.embed_dir, 'altsplice')
    EVENT_LIST = ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime', 'concatenated']

    TSNE_PP_PLOT_SET  = [50]
    TSNE_LR_PLOT_SET  = [500]
    do_pdf = True

    nan_class = 3
    if np.any(gtdf.isnull()):
        assert not nan_class in gtdf.values
        gtdf.values[np.isnan(gtdf.values)] = nan_class
    gtdf = gtdf.apply(lambda x: pd.to_numeric(x, downcast='integer'))
    mut_cmap = sns.color_palette('Set2', nan_class+1)

    for event in EVENT_LIST:
        embed_base_dir = os.path.join(config.embed_dir, 'altsplice', event)
        pca_model, pca_embeds, tsne_embeds_dict = utils.load_embeds(embed_base_dir)
        pca_embeds.index = pca_embeds.index.map(lambda x: x.replace('.aligned', ''))
        pca_embeds.index = pca_embeds.index.map(lambda x: x.replace('.npz', ''))
        for df in tsne_embeds_dict.values():
            df.index = df.index.map(lambda x: x.replace('.aligned', ''))
            df.index = df.index.map(lambda x: x.replace('.npz', ''))
        metadata_df = utils.load_metadata_df(config.metadata_path, pca_embeds.index)

        # Define plotting kwargs
        default_kwargs = {'alpha': .25, 's':10}
        plotting_df = pd.DataFrame(index=metadata_df.index)
        assert isinstance(mut_cmap[nan_class], tuple)
        plotting_df['facecolors'] = pd.Series([mut_cmap[nan_class]], index=plotting_df.index)
        plotting_df['marker'] = pd.Series('o', index=plotting_df.index)

        # List of kwargs used to create a legend
        legend_kwargs_list = list()
        legend_defaults = {'lw':0, 's': 30, 'alpha':.75}
        for i, color in enumerate(mut_cmap):
            label = str(i) if i != nan_class else 'Missing'
            legend_kwargs_list.append(_copy_and_update(legend_defaults, {'c':color, 'label': label}))

        embed_base_out = os.path.join(config.plot_dir, 'altsplice', event, 'embeds', 'sqtl')
        for gene in gtdf:
            ensg = gene.split('_')[0].split('.')[0]
            mut_status = gtdf[gene]
            tag = gene.replace('.','-')
            for (pp, lr), tsne_embeds in tsne_embeds_dict.items():
                if not TSNE_PP_PLOT_SET == 'all' and not pp in TSNE_PP_PLOT_SET: continue
                if not TSNE_LR_PLOT_SET == 'all' and not lr in TSNE_LR_PLOT_SET: continue
                tsne_plot_out = os.path.join(embed_base_out, ensg, '%s_tsne_embeds_pp_%d_lr_%d.png' %(tag, pp, lr))
                if not os.path.exists(os.path.dirname(tsne_plot_out)): os.makedirs(os.path.dirname(tsne_plot_out))
                plotting_df['facecolors'] = pd.Series([mut_cmap[nan_class]], index=plotting_df.index)
                plotting_df['facecolors'].loc[mut_status.index] = pd.Series(mut_status.map(lambda x: mut_cmap[x]), index=mut_status.index)
                fig, ax = ebhelp.plot_tsne_embeddings(tsne_embeds, plotting_df, legend_kwargs_list, default_kwargs)
                fig.suptitle('AltSplice %s tSNE(perplexity=%d, learning_rate=%d)' %(tag, pp, lr))
                print "Writing %s" %tsne_plot_out
                plt.savefig(tsne_plot_out, bbox_inches='tight', dpi=300)
                if do_pdf: ebhelp.save_as_pdf(tsne_plot_out)
                plt.close()
    return

MAX_EVENTS = 5000
DEBUG = True

if __name__ == '__main__':
    gtdf = load_gt_data()

    # event_type to path of raw psi values
    gtdf = gtdf.apply(lambda x: pd.to_numeric(x, downcast='integer'))
    gtdf = gtdf.apply(lambda x: pd.to_numeric(x, downcast='integer'))
    map_etype_to_file = {'exon_skip': config.alt_splice_exon_skip_path,
                         'intron_retention': config.alt_splice_intron_retention_path,
                         'alt_3prime': config.alt_splce_alt_3prime_path,
                         'alt_5prime': config.alt_splce_alt_5prime_path}
    etype_desc = 'combined'
#    do_heatmap(combined_df, gtdf, col_cmap_lut, combined_df_desc)
#    do_embeddings(gtdf)

