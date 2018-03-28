import os
import sys
import re
import numpy as np
import scipy as sp
import pandas as pd
import h5py
import fnmatch
import os

sys.path.append('/cluster/home/starks/git/tools-python/viz')
import gwas as plotter

BASEDIR = os.path.dirname(os.path.dirname(__file__))
sys.path.append(BASEDIR)
import config
import utils

def _load_bookmark(reset):
    path = os.path.join(BKDIR, 'bookmark')
    if reset or not os.path.exists(path):
        return 0
    assert THOLD == .01 / 300000
    return int(open(path, 'r').read().strip())

def _write_bookmark(bookmark):
    print "Saving progress, bookmark = %d" %bookmark
    path = os.path.join(BKDIR, 'bookmark')
    open(path, 'w').write(str(bookmark))
    return


def _load_skip_paths(reset):
    assert THOLD == .01 / 300000
    skip_paths = set()
    path = os.path.join(BKDIR, 'below_thold.paths')
    if reset:
        open(path, 'w').close()
    if os.path.exists(path):
        skip_paths.update(open(path, 'r').read().splitlines())
    return skip_paths

def _write_skip_paths(skip_paths):
    path = os.path.join(BKDIR, 'below_thold.paths')
    with open(path, 'w') as out:
        out.write('\n'.join(list(skip_paths)))
    return

def name_figpath(hdf5_path, gene_name=None):
    outpath = os.path.join(PLOT_DIR, os.path.relpath(hdf5_path, WDIR))
    outpath = re.sub('.hdf5$', '.png', outpath)
    if not gene_name is None:
        dirname = os.path.dirname(outpath)
        basename = os.path.basename(outpath)
        outpath = os.path.join(dirname, 'gene_%s_%s'%(gene_name, basename))
    return outpath

def save_as_pdf(png_path):
    assert png_path.endswith('.png')
    pdf_path = re.sub('.png$', '.pdf', png_path)
    print "Writing: %s"%pdf_path
    plotter.plt.savefig(pdf_path)
    return

def process(path, gene_name):
    outpath = name_figpath(path, gene_name)
    fin = h5py.File(path, 'r')
    pv = fin['pv'][:]
    if sp.nanmin(pv) < THOLD:
        pos = fin['pos'][:]
        plotter.makeManhattanPlot(pv, pos, thold_denom=THOLD_DENOM)
        if not os.path.exists(os.path.dirname(outpath)): os.makedirs(os.path.dirname(outpath))
        print "Writing %s" %outpath
        plotter.plt.savefig(outpath, dpi=100)
        if WRITE_PDF_PLOTS: save_as_pdf(outpath)
        plotter.plt.close()
    return

def _load_census_genes(path):
    '''Loads genes (ENSGXXXXX... format).
    '''
    genes = np.loadtxt(path, dtype=str)
    return genes

def _create_single_map_gene_idx_to_census_gene(gene_idx, gene_names, census_genes):
    '''Creates a map from gene_idx to census genes.

    gene_idx: maps genes in gene_names to idx WARNING: ASSUMES INDEX STARTS AT 1
    gene_names: genes in ENSG form, may have .# ending
    census_genes: genes in ENSG form, .# ending clipped

    Returns
    gene_idx_mask: bool mask over gene_idx for inclusion in census genes
    masked_gene_names: name for each idx in mask
    '''
    assert gene_idx.max() < gene_names.size
    clipped_gene_names = map(lambda x: x.split('.')[0], gene_names)
    censor_idxs = np.where(np.in1d(clipped_gene_names, census_genes))[0]
    gene_idx_mask = np.in1d(gene_idx-1, censor_idxs) # HERE IS THE FIX FOR INDEX STARTING AT 1
    masked_gene_names = gene_names[gene_idx[gene_idx_mask]]
    assert masked_gene_names.size == gene_idx_mask.sum()
    gene_idx_in_census = np.where(gene_idx_mask)[0]
    map_gene_idx_to_census_gene = dict(zip(gene_idx_in_census, masked_gene_names))
    assert len(map_gene_idx_to_census_gene) > 0
    return map_gene_idx_to_census_gene

def _create_maps_gene_idx_to_census_genes(census_genes, pheno_h5_fin):
    '''Matches gene_idx rows to census_genes of interest.

    Arguments
    census_genes: array of gene names (ENSG format)
    pheno_h5_fin: hdf5 file with keys: [alt_3prime, alt_5prime, etc]

    Returns:
    A dict with same keys as pheno_h5_fin, value is a mask over gene_idx.
    '''
    dict_map_gene_idx_to_census_gene = dict()
    for key in pheno_h5_fin.keys():
        if not key in _EVENT_TYPES: continue
        data = pheno_h5_fin[key]
        map_gene_idx_to_census_gene = _create_single_map_gene_idx_to_census_gene(
                                          data['gene_idx'][:],
                                          data['gene_names'][:],
                                          census_genes)
        dict_map_gene_idx_to_census_gene[key] = map_gene_idx_to_census_gene
    return dict_map_gene_idx_to_census_gene

def _extract_event_path_fields(event_path):
    bname = os.path.basename(event_path)
    fields = bname.split('_')
    assert len(fields) == 12
    assert fields[0] == 'chunk'
    assert fields[-1] == 'germ.hdf5'
    idx = int(fields[1])
    event_pos = tuple(map(lambda x: int(x), fields[5:-1]))
    assert len(event_pos) == 6
    return idx, event_pos

def _create_map_path_to_etype(list_of_event_paths):
    '''Extracts the event from the path.

    Assumes paths are in format: $WDIR/$event/...
    '''
    etypes = map(lambda x: os.path.relpath(x, WDIR).split('/')[0], list_of_event_paths)
    assert np.all(np.in1d(np.unique(etypes), _EVENT_TYPES))
    return dict(zip(list_of_event_paths, etypes))

def _event_pos_is_matched(idx, event_pos, etype_pheno_h5_fin):
    return np.array_equal(event_pos, etype_pheno_h5_fin['event_pos'][:, idx])

def _data_exists_on_all_events(map_path_to_etype, splicing_pheno_h5_fin):
    unique_etypes = np.unique(map_path_to_etype.values())
    return np.all(np.in1d(unique_etypes, splicing_pheno_h5_fin.keys()))

def _do_check_idx_matches():
    return np.random.rand() < .05

def filter_to_census_genes(list_of_event_paths, census_genes, splicing_pheno_h5_fin):
    '''Filter paths to events that occur within census genes.

    list_of_event_paths: list of all hdf5 paths, (holding sqtl data)
    census_genes: list of genes in the census (ENSG format)
    splicing_pheno_h5_fin: holds data used to generate event hdf5 paths
    '''
    all_etypes = [u'alt_3prime', u'alt_5prime', u'exon_skip', u'intron_retention']
    map_path_to_etype = _create_map_path_to_etype(list_of_event_paths)
    assert _data_exists_on_all_events(map_path_to_etype, splicing_pheno_h5_fin)
    dict_map_gene_idx_to_census_genes = _create_maps_gene_idx_to_census_genes(census_genes, splicing_pheno_h5_fin)

    census_paths = list()
    map_path_to_gene_name = dict()
    for path, etype in map_path_to_etype.iteritems():
        etype_pheno_h5_fin = splicing_pheno_h5_fin[etype]
        idx, event_pos = _extract_event_path_fields(path)
        if _do_check_idx_matches():  assert _event_pos_is_matched(idx, event_pos, etype_pheno_h5_fin)
        gene_name = dict_map_gene_idx_to_census_genes[etype].get(idx)
        if gene_name is None: continue

        map_path_to_gene_name[path] = gene_name

    return map_path_to_gene_name

THOLD_DENOM = 300000
THOLD = .05 / THOLD_DENOM

_EVENT_TYPES = [u'alt_3prime', u'alt_5prime', u'exon_skip', u'intron_retention']
WDIR = '/cluster/work/grlab/projects/TCGA/PanCancer/QTL_Analysis/transAssociation_10K/sQTL/sc_pancanatlas_r10_s2000_V100_old'
BKDIR = os.path.join(WDIR, '.bookkeep')
#SKIP_PATHS = _load_skip_paths(True)
PLOT_DIR = os.path.join(config.plot_dir, 'altsplice', 'sqtl_manhattan')
WRITE_PDF_PLOTS = False
OVERWRITE = True
#BOOKMARK = _load_bookmark(reset=False)
DEBUG = True

if __name__ == '__main__':
    event_h5_path_list = open(os.path.join(BKDIR, 'hdf5.paths'), 'r').read().splitlines()
    census_genes = _load_census_genes(config.census_gene_path)
    splicing_pheno_h5_path = '/cluster/work/grlab/projects/TCGA/PanCancer/hdf5_10Kset/phenotypes.hdf5r10_s2000_V100.hdf5'
    splicing_pheno_h5_fin = h5py.File(splicing_pheno_h5_path, 'r')['splicing']
    map_event_h5_paths_to_census_gnames = filter_to_census_genes(event_h5_path_list, census_genes, splicing_pheno_h5_fin)
    for path, gname in map_event_h5_paths_to_census_gnames.iteritems():
        process(path, gname)

