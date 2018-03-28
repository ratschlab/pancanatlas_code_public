import os
import re
import h5py
import pandas as pd
import numpy as np
import glob
import re
import config

_GREY_HEX = '#D3D3D3'

def parse_tcga_id(tcga_id):
    tcga, tss, participant, sample, portion, plate, center = tcga_id.split('-')[:7]
    assert tcga.lower() == 'tcga'
    is_tumor = int(sample[0]) == 0
    return tss, is_tumor

def get_tss(tcga_id):
    return tcga_id.split('-')[1]

def tcga_id_is_tumor(tcga_id):
    assert tcga_id.upper().startswith('TCGA')
    tumor_tag = int(tcga_id.split('-')[3][0])
    is_tumor = tumor_tag == 0
    return is_tumor

def load_embeds(embed_dir, whitelist=None, pp_set='all', lr_set='all'):
    '''Returns pca_model, pca_embeds, tsne_embeds_dict
    '''
    # Load embeddings, models
    pca_embeds_path_glob = glob.glob(os.path.join(embed_dir, 'pca_embed*.tsv'))
    assert len(pca_embeds_path_glob) == 1
    pca_embeds_path = pca_embeds_path_glob[0]
    pca_model_base = os.path.basename(pca_embeds_path).replace('embeds', 'model').replace('.tsv', '.npy')
    pca_model_path = os.path.join(embed_dir, pca_model_base)


    print "Loading PCA embeddings from %s" %pca_embeds_path
    pca_embeds = pd.read_csv(pca_embeds_path, sep='\t', index_col=0)
    pca_model = np.load(pca_model_path).item()
    pca_embeds.index = [re.sub('.aligned$', '', x) for x in pca_embeds.index]
    if whitelist is not None: pca_embeds = pca_embeds.loc[whitelist].dropna()

    tsne_embeds_dict = dict()
    _pp_searcher = re.compile('pp_([0-9]+)')
    _lr_searcher = re.compile('lr_([0-9]+)')
    for tsne_path in glob.glob(os.path.join(embed_dir, 'tsne_embeds*.tsv')):
        pp = int(_pp_searcher.search(tsne_path).group(1))
        lr = int(_lr_searcher.search(tsne_path).group(1))
        if not pp_set == 'all' and not pp in pp_set: continue
        if not lr_set == 'all' and not lr in lr_set: continue
        key = (pp, lr)
        tsne_embeds = pd.read_csv(tsne_path, sep='\t', index_col=0)
        tsne_embeds.index = [re.sub('.aligned$', '', x) for x in tsne_embeds.index]
        if whitelist is not None: tsne_embeds = tsne_embeds.loc[whitelist].dropna()
        tsne_embeds_dict[key] = tsne_embeds
    return pca_model, pca_embeds, tsne_embeds_dict

def load_libsize(libsize_path):
    return pd.read_csv(libsize_path, sep='\t', index_col=0)

def load_rnadeg(rnadeg_path):
    df = pd.read_csv(rnadeg_path, sep='\t', index_col=0, header=None, names=['RNADeg_Score'])
    df.index = [os.path.basename(pth.replace('.aligned.hdf5', '')) for pth in df.index]
    return df

def load_whitelist(whitelist_path):
    print("Loading whitelist from %s" %whitelist_path)
    return np.loadtxt(whitelist_path, dtype='str')

def load_cnc_series(config):
    # Loading metadata
    print "Loading metadata from %s" %config.processed_expression_count_path
    proc_expression_file = h5py.File(config.processed_expression_count_path, 'r')
    cnc = proc_expression_file['cnc'][:]
    # plate = proc_expression_file['plate'][:]
    gtids = proc_expression_file['gtids'][:]
    cnc_series = pd.Series(cnc, index=gtids)
    cnc_series = cnc_series.loc[cnc_series != 'NA'].drop('NA', errors='ignore')
    proc_expression_file.close()
    return cnc_series

def load_cnc_and_is_tumor_series(config):
    cnc_series = load_cnc_series(config)
    is_tumor_series = pd.Series([tcga_id_is_tumor(tcga_id) for tcga_id in cnc_series.index], index=cnc_series.index)
    return cnc_series, is_tumor_series


def load_metadata_df(path, embeds_index=None):
    '''Loads metadata, constrains to ids in embeds_index if given.
    '''
    df = pd.read_csv(path, sep='\t')[['tcga_id', 'study', 'is_normal']]
    df = df.drop_duplicates('tcga_id').set_index('tcga_id')
    df.iloc[:]['is_normal'] = ~df.iloc[:]['is_normal'].values
    df.columns = ['cnc', 'is_tumor']
    if not embeds_index is None:
        df = translate_tcga_to_strain_index(df, embeds_index)
        df = df.loc[df.index & embeds_index]
    return df

def translate_tcga_to_strain_index(metadata_df, index):
    metadata_df_trans = pd.DataFrame(index=index, columns=metadata_df.columns)
    for x in index:
       metadata_df_trans.loc[x] = metadata_df.loc[x.split('.')[0]]
    for col, dtype in metadata_df.dtypes.to_dict().items():
        metadata_df_trans[col] = metadata_df_trans[col].astype(dtype)
    return metadata_df_trans

def load_color_scheme(path):
    df = pd.read_csv(path, sep='\t')[['Hex Colors', 'Study Abbreviation']]
    cnc_to_color = dict(zip(df['Study Abbreviation'].values, df['Hex Colors'].values))
    return cnc_to_color

def append_subtype(path, metadata_df):
    # do not match normals
    df = pd.read_csv(path, sep='\t').set_index('pan.samplesID')
    index = metadata_df.loc[metadata_df['is_tumor']].index
    map_index_old_new = dict()
    for dfx in df.index:
        matches = map(lambda x: x.startswith(dfx), index.values)
        if sum(matches) == 1: map_index_old_new[dfx] = index[matches][0]
    #    if sum(matches) > 1: print dfx
    df = df.loc[map_index_old_new.keys(), [st for st in df.columns if st.startswith('Subtype')]]
    subtype_names = df.columns
    df.index = df.index.map(lambda x: map_index_old_new[x])
    df = pd.concat([metadata_df, df], axis=1)
    return df, subtype_names

def load_staging(path, index):
    '''Load staging data, matching to index
    '''
    df = pd.read_csv(path, sep='\t')
    df = df.replace(re.compile('\[*\]'), np.nan).set_index('bcr_patient_barcode')
    df = df.drop(['bcr_patient_uuid', 'acronym'], axis=1)
    index_map = pd.Series(index.values, index=index.map(lambda x: '-'.join(x.split('-')[:3])))
    filtered_index = index_map.loc[df.index].dropna()
    df = df.loc[filtered_index.index]
    df.index = filtered_index.values
    return df

def append_staging(path, metadata_df):
    staging_df = load_staging(config.staging_info, metadata_df.index)
    subtype_names = staging_df.columns
    metadata_df = metadata_df.loc[staging_df.index]
    df = pd.concat((metadata_df, staging_df), axis=1)
    return df, subtype_names

def load_smoking_and_age(path, index):
    df = pd.read_csv(path, sep='\t')
    df = df.replace(re.compile('\[*\]'), np.nan).set_index('bcr_patient_barcode')
    df = df.drop(['bcr_patient_uuid', 'acronym'], axis=1)
    index_map = pd.Series(index.values, index=index.map(lambda x: '-'.join(x.split('-')[:3])))
    filtered_index = index_map.loc[df.index].dropna()
    df = df.loc[filtered_index.index]
    df.index = filtered_index.values
    return df

def append_smoking_and_age(path, metadata_df):
    staging_df = load_staging(path, metadata_df.index)
    subtype_names = staging_df.columns
    metadata_df = metadata_df.loc[staging_df.index]
    df = pd.concat((metadata_df, staging_df), axis=1)
    return df, subtype_names

def get_cache_base(path):
    path = os.path.abspath(path)
    assert path.startswith(config.PROJECT_DIR)
    cache = re.sub('^%s/'%config.PROJECT_DIR, '', path)
    return os.path.join(config.CACHE_DIR, os.path.splitext(cache)[0])

def cache_large_df(basename, df):
    index_path = basename + '_index.npy'
    np.save(index_path, df.index)
    values_path = basename + '_values.npy'
    np.save(values_path, df.values)
    columns_path = basename + '_columns.npy'
    np.save(columns_path, df.columns)
    return

def load_large_df(basename):
    index_path = basename + '_index.npy'
    index = np.load(index_path)
    values_path = basename + '_values.npy'
    values = np.load(values_path)
    index = np.load(index_path)
    columns_path = basename + '_columns.npy'
    columns = np.load(columns_path)
    return pd.DataFrame(values, index=index, columns=columns)



