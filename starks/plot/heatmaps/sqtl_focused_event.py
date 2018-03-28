
_HEATMAP_BASE = os.path.join(config.plot_dir, 'altsplice', 'concatenated', 'heatmaps', 'sqtl')
_EMBED_BASE = os.path.join(config.plot_dir, 'altsplice', 'concatenated', 'embeds', 'sqtl')

_TRANS_DATA_PATH = '/cluster/work/grlab/projects/TCGA/PanCancer/QTL_Analysis_results/analysis_10kset/2dmanhattan_information_v3/transdataresults_test_v3.tsv'
_GENOTYPE_H5_PATH  = '/cluster/work/grlab/projects/TCGA/PanCancer/hdf5_10Kset/mergedFiles_clean.hdf5'


def load_trans_info(trans_path):
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

def find_genes_with_targets(trans_df, k=50):
    targets = trans_df.loc[:, 'gene-id(mutation)']
    genes, counts = np.unique(targets.values, return_counts=True)
    return genes[counts >= k]

def load_mutation_info(trans_df, genes_with_targets):
    # get all coords
    # for all coords get mutation information for all samples
    coords = get_coords(trans_df, genes_with_targets)
    gt_file = h5py.File(sqtl._GENOTYPE_H5_PATH, 'r')
    gtdf = sqtl.get_mutation_info(gt_file, coords)
    return gtdf


def get_coords(trans_df, target_genes):
    target_mask = trans_df.loc[:, 'gene-id(mutation)'].isin(target_genes)
    coords = trans_df.loc[target_mask, ['snp_chrm', 'snp_pos', 'gene-id(mutation)']]
    coords = coords.drop_duplicates()
    coords.index = np.arange(coords.shape[0])
    coords.columns = ['chr', 'pos', 'gene']
    return coords


