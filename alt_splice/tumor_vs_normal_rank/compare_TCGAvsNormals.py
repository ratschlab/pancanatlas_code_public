import sys
import os
import pdb
import re
import gzip
import scipy as sp
import scipy.stats as spst
import scipy.io as scio
import h5py
import cPickle
import glob
import numpy.random as npr
npr.seed(23)
import fisher

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pylab import cm
import matplotlib.gridspec as gridspec

sys.path.append('/cluster/home/akahles/git/tools/python/viz')
from distribution import violin_plot
from distribution import dist_overview
sys.path.append('/cluster/home/akahles/git/tools/python/oop_libs')
from expression import ExpressionData

sys.path.append('/cluster/home/akahles/git/projects/2013/PanCancerTCGA/rerun2017')
from utils.paths import paths
import utils.utils as utils
import utils.samples as samples
import utils.names as names

###
### Directories and data access
###
basedir_anno = '/cluster/work/grlab/projects/TCGA/PanCancer/rerun_alt_splice_anno'
metadata_style = 'pancan_rerun18'

###
### Settings
###
event_types = ['exon_skip', 'alt_3prime']
CONF = 2

class options(object):
    """Class that contains all the optional parameters
       to be passed around."""
    pass

if len(sys.argv) < 2:
    options.cancer_type = 'all'
else:
    options.cancer_type = sys.argv[1]

# path to file that contains sample IDs to be excluded from analysis
# file should contain one ID per line
# e.g. options.blacklist = '/cluster/project/grlab/home/akahles/git/projects/2014/icgc/alternative_splicing/sample_blacklist'
options.blacklist = None        # path to file that 

# list of genes to be excluded from analysis
# e.g. options.gene_blacklist = ['ENSG00000164304.11', 'ENSG00000073754.5']
options.gene_blacklist = []

# max avg fraction of samples that have evidence
# dictionary with study ID as key and float (between 0 and 1) as value
options.max_norm_expression_frac = {'tn':0.01,
                                    'gt':0.01}

# restrict analysis to a given set of cancer census genes
options.census_only = False
options.census_file = '../annotation/census'

# lower count thresholds to recognize an intron as expressed in a sample
# dictionary with study ID as key and min read count (int) as value
options.thresholds = {'tc':10,
                      'tn':3,
                      'gt':2}

# remove introns that are higher expressed in normals than in tumors
# (computed individually over the samples showing the intron)
options.filter_high_normal_expression = True

### minimal number of expressed introns for a sample to be included
options.min_introns_per_sample = 5

# number of top events to be considered for the plot
options.top = 40

# definition of tumor and normal sample sets
options.tumor_set = ['tc']
options.normal_set = ['tn', 'gt']#, 'ec']

# in which sets should the IDs be shortened? Only relevant for TCGA sample IDs.
options.shorten_label = [('tc', 1), ('tn', 1), ('gt', 0)]
#options.shorten_label = []

options.plot_intron_distribution = False

# cancer types to be considered for analysis
cancer_type_tag = ''
if options.cancer_type == 'all':
    which_types = {'tc':['BRCA', 'COAD', 'LIHC', 'LUAD', 'LUSC', 'PRAD', 'THCA'], ### cancer types with at least 40 normals and 40 GTEx samples
                   'tn':['BRCA', 'COAD', 'LIHC', 'LUAD', 'LUSC', 'PRAD', 'THCA'],
                   'gt':['BRST', 'COLN', 'LIVR', 'LUNG', 'PRST', 'THYR']} #, 'XTR.']}
    options.width = 3
elif options.cancer_type == 'BLCA':
    which_types = {'tc':['BLCA'], 'tn':['BLCA'], 'gt':['BLDR']}
    cancer_type_tag = '_BLCA'
    options.width = 2
elif options.cancer_type == 'BRCA':
    which_types = {'tc':['BRCA'], 'tn':['BRCA'], 'gt':['BRST']}
    cancer_type_tag = '_BRCA'
    options.width = 2
elif options.cancer_type == 'CHOL':
    which_types = {'tc':['CHOL'], 'tn':['CHOL'], 'gt':['LIVR']}
    cancer_type_tag = '_CHOL'
    options.width = 2
elif options.cancer_type == 'COAD':
    which_types = {'tc':['COAD'], 'tn':['COAD'], 'gt':['COLN']}
    cancer_type_tag = '_COAD'
    options.width = 2
elif options.cancer_type == 'ESCA':
    which_types = {'tc':['ESCA'], 'tn':['ESCA'], 'gt':['STOM']}
    cancer_type_tag = '_ESCA'
    options.width = 2
elif options.cancer_type == 'READ':
    which_types = {'tc':['READ'], 'tn':['READ'], 'gt':['COLN']}
    cancer_type_tag = '_READ'
    options.width = 2
elif options.cancer_type == 'HNSC':
    which_types = {'tc':['HNSC'], 'tn':['HNSC']}
    cancer_type_tag = '_HNSC'
    options.width = 2
elif options.cancer_type == 'KICH':
    which_types = {'tc':['KICH'], 'tn':['KICH'], 'gt':['KIDN']}
    cancer_type_tag = '_KICH'
    options.width = 2
elif options.cancer_type == 'KIRC':
    which_types = {'tc':['KIRC'], 'tn':['KIRC'], 'gt':['KIDN']}
    cancer_type_tag = '_KIRC'
    options.width = 2
elif options.cancer_type == 'KIRP':
    which_types = {'tc':['KIRP'], 'tn':['KIRP'], 'gt':['KIDN']}
    cancer_type_tag = '_KIRP'
    options.width = 2
elif options.cancer_type == 'LIHC':
    which_types = {'tc':['LIHC'], 'tn':['LIHC'], 'gt':['LIVR']}
    cancer_type_tag = '_LIHC'
    options.width = 2
elif options.cancer_type == 'LUAD':
    which_types = {'tc':['LUAD'], 'tn':['LUAD'], 'gt':['LUNG']}
    cancer_type_tag = '_LUAD'
    options.width = 2
elif options.cancer_type == 'LUSC':
    which_types = {'tc':['LUSC'], 'tn':['LUSC'], 'gt':['LUNG']}
    cancer_type_tag = '_LUSD'
    options.width = 2
elif options.cancer_type == 'PRAD':
    which_types = {'tc':['PRAD'], 'tn':['PRAD'], 'gt':['PRST']}
    cancer_type_tag = '_PRAD'
    options.width = 2
elif options.cancer_type == 'STAD':
    which_types = {'tc':['STAD'], 'tn':['STAD'], 'gt':['STOM']}
    cancer_type_tag = '_STAD'
    options.width = 2
elif options.cancer_type == 'THCA':
    which_types = {'tc':['THCA'], 'tn':['THCA'], 'gt':['THYR']}
    cancer_type_tag = '_THCA'
    options.width = 2
elif options.cancer_type == 'UCEC':
    which_types = {'tc':['UCEC'], 'tn':['UCEC'], 'gt':['UTER']}
    cancer_type_tag = '_UCEC'
    options.width = 2

plot_tag = ''

options.use_whitelist = True
wl_tag = ''
if options.use_whitelist:
    wl_tag = '.whitelisted'
plot_tag += wl_tag

options.subsample = True # False
options.subsample_size = 35
subs_tag = ''
if options.subsample:
    subs_tag = '.subsampled%i' % options.subsample_size
plot_tag += subs_tag

### when testing per tumor type we will use the following dictionary
### to assign tumor normal pairs
test_dict = dict()
test_dict['gt'] = {'BRCA-US':'BRST',
                   'KICH-US':'KIDN',
                   'KIRC-US':'KIDN',
                   'LIHC-US':'LIVR',
                   'LIRI-JP':'LIRI',
                   'LUAD-US':'LUNG',
                   'RECA-EU':'KIDN'}

# if true, the introns are extracted from the graph rather than from 
# the individual events
graph_introns = True

# if true, all introns that also ocurr in the annotation are excluded
no_anno = False

# experimental - if true, restrict to events/genes that also are sQTL
options.sqtl_only = False
options.sqtl_file = ''

# experimental - plot the expression distribution? TODO!
plot_distributions = False

# choose the strategy to rank the tumor specific events
#   sum_non_tumor:
#   fisher_subset:
#   ...
options.strategy = 'fisher_subset'
fisher_subset = {'tn':['BRCA', 'COAD', 'LIHC', 'LUAD', 'LUSC', 'PRAD', 'THCA'],
                 'gt':['BRST', 'COLN', 'LIVR', 'LUNG', 'PRST', 'THYR']}

# output directory for plots
options.plotdir = os.path.join(paths.plotdir, 'tss_introns', 'plots_tumor_specific_splicing%s_%s' % (cancer_type_tag, options.strategy))
options.plotdir = utils.gen_plot_directory(options)
if not os.path.exists(options.plotdir):
    os.makedirs(options.plotdir)

###
### Peliminaries
###

# get lookup table of gene names
lookup = names.get_lookup_complete()

# create cancer type dictionary
(ct_dict, is_tumor_dict) = utils.get_ct_dict_metatable(paths.metadata, style=metadata_style)

# create GTex type dictionary
sample_dict = {'BLDR':'../annotation/sample_lists_gtex/Bladder.txt',
               'BRST':'../annotation/sample_lists_gtex/Breast.txt',
               'COLN':'../annotation/sample_lists_gtex/Colon.txt',
               'KIDN':'../annotation/sample_lists_gtex/Kidney.txt',
               'LIVR':'../annotation/sample_lists_gtex/Liver.txt',
               'LUNG':'../annotation/sample_lists_gtex/Lung.txt',
               'PRST':'../annotation/sample_lists_gtex/Prostate.txt',
               'STOM':'../annotation/sample_lists_gtex/Stomach.txt',
               'THYR':'../annotation/sample_lists_gtex/Thyroid.txt',
               'UTER':'../annotation/sample_lists_gtex/Uterus.txt',
               'CELL':'../annotation/sample_lists_gtex/Cells.txt'}
gt_dict = utils.get_gt_dict(sample_dict)


###
### Actual work
###

### load expression sample from hdf5 files
exp_hdf5s = [os.path.join(paths.basedir, 'rerun2018_hdf5', 'expression_counts%s.hdf5' % wl_tag),
             os.path.join(paths.basedir_gtex, 'counts_gtex_2015-03-29.hdf5')]

exp_tumor_rule = [re.compile(r'TCGA-..-....-0.*'),
                  re.compile(r'XXXXXX')]
exp_sample = samples.expression_sample.from_hdf5(exp_hdf5s, rules=exp_tumor_rule)
exp_sample.assign_tumor_samples(is_tumor_dict)

tumor_sids = exp_sample.sids[sp.where(exp_sample.is_tumor)[0]]
normal_sids = exp_sample.sids[sp.where(~exp_sample.is_tumor)[0]] 

### load intron counts from HDF5
### list of hdf5 files to take into account
hdf5_files = {'tc' : os.path.join(paths.basedir_as, 'spladder', 'genes_graph_conf%i.merge_graphs.validated.count%s.hdf5' % (CONF, wl_tag)),
              'tn' : os.path.join(paths.basedir_as, 'spladder', 'genes_graph_conf%i.merge_graphs.validated.count%s.hdf5' % (CONF, wl_tag)),
              'gt' : os.path.join(paths.basedir_as_gtex, 'spladder', 'genes_graph_conf%i.merge_graphs.validated.count.hdf5' % CONF)}
              
hdf5_rules = {'tc' : tumor_sids,
              'tn' : normal_sids}

anno = os.path.join(paths.basedir_as, 'spladder', 'genes_graph_conf%i.merge_graphs.validated.pickle' % CONF)
            
intron_samples = utils.preprocess_intron_counts(None, paths.basedir_as, hdf5_files, hdf5_rules, annotation=anno, dataset='pancanatlas', tag=wl_tag)

print >> sys.stdout, 'Found %i unique introns' % intron_samples[options.tumor_set[0]].eid.shape[0]

### blacklist samples as necessary
if options.blacklist is not None:
    bl_data  = sp.loadtxt(options.blacklist, dtype='str', delimiter='\t')#[:, 0]
    for t in intron_samples:
        intron_samples[t].subset_by_strain(bl_data, inverse=True)
    exp_sample.subset_by_strain(bl_data, inverse=True)

### shorten ids if necessary
for t in options.shorten_label:
    intron_samples[t[0]].strains = sp.array([x.split('.')[t[1]] for x in intron_samples[t[0]].strains])
exp_sample.sids = sp.array([x.split('.')[0] if x.startswith('SRR') else x.split('.')[1] for x in exp_sample.sids], dtype='str')

### make unique over samples
for t in intron_samples:
    intron_samples[t].make_unique_by_strain()

### sort TCGA data by cancer/tissue type
intron_samples['tc'].add_types(ct_dict) ## replaces ctypes_t
intron_samples['tc'].sort_by_type()
if 'tn' in intron_samples:
    intron_samples['tn'].add_types(ct_dict) ## replaced ctypes_n
    intron_samples['tn'].sort_by_type()
if 'gt' in intron_samples:
    ### sort GTex samples by tissue type
    intron_samples['gt'].add_types(gt_dict)
    intron_samples['gt'].sort_by_type()
    ### remove samples that could not be matched to tissue type
    intron_samples['gt'].remove_unmatched_types()

### subset to cancer types where we have normal data
if 'tn' in intron_samples:
    k_idx = sp.where(sp.in1d(intron_samples['tc'].types, intron_samples['tn'].types))[0]
    intron_samples['tc'].reorder_by_col_idx(k_idx)

### subset to a given set of cancer/tissue types
for t in which_types:
    intron_samples[t].subset_by_type(which_types[t])

### perform subsampling
if options.subsample:
    for t in intron_samples:
        k_idx = []
        for tt in sp.unique(intron_samples[t].types):
            tt_idx = sp.where(intron_samples[t].types == tt)[0]
            if tt_idx.shape[0] > options.subsample_size:
               tt_idx = sp.sort(npr.choice(tt_idx, size=options.subsample_size, replace=False))
            k_idx.extend(tt_idx)
        k_idx = sp.sort(sp.array(k_idx))
        intron_samples[t].reorder_by_col_idx(k_idx)

### remove introns that are not expressed in tumor
k_idx = sp.where(sp.sum(intron_samples['tc'].introns, axis=1) > 0)[0]
print >> sys.stdout, 'Removed %i introns that showed no expression in tumor set' % (intron_samples['tc'].eid.shape[0] - k_idx.shape[0])
print >> sys.stdout, 'Retaining %i introns' % k_idx.shape[0]
for t in intron_samples:
    intron_samples[t].reorder_by_row_idx(k_idx)

### remove by gene blacklist
if len(options.gene_blacklist) > 0:
    intron_samples.subset_by_gene_id(options.gene_blacklist, inverse=True)

### subset to introns in associated genes only
sqtl_tag = ''
if options.sqtl_only:
    sqtl_tag = '.sqtl'
    assoc = sp.loadtxt(options.sqtl_file, dtype='str')
    gid_short = sp.array([x.split('.')[0] for x in intron_samples['tc'].gid])
    k_idx = sp.where(sp.in1d(gid_short, assoc))[0]
    for t in intron_samples:
        intron_samples[t].reorder_by_row_idx(k_idx)
plot_tag += sqtl_tag

### remove annotated introns
anno_tag = ''
if no_anno:
    ### generate anno event ID strings
    print 'Generating anno event IDs ...'
    anno_genes = cPickle.load(open(os.path.join(basedir_anno, 'spladder', 'genes_graph_conf3.merge_graphs.pickle'), 'r'))[0]
    ### reduce events to inner coordinates 
    anno_ids = []
    for x in anno_genes:
        for tr in x.exons:
            if tr.shape[0] < 2:
                continue
            for e in range(1, tr.shape[0]):
                if tr[e-1, 1] <  tr[e, 0]:
                    anno_ids.append('.'.join([x.chr, str(tr[e-1, 1]), str(tr[e, 0])]))
                else:
                    anno_ids.append('.'.join([x.chr, str(tr[e, 1]), str(tr[e-1, 0])]))
    anno_ids = sp.unique(sp.array(anno_ids))

    ### subset intron list to unannotated ones
    for t in intron_samples:
        intron_samples[t].subset_by_event_id(anno_ids, inverse=True)
    anno_tag = '.no_anno'
plot_tag += anno_tag

### remove introns that are not in cancer census genes
census_tag = ''
if options.census_only:
    census = sp.loadtxt(options.census_file, dtype='str')
    for t in intron_samples:
        intron_samples[t].subset_by_gene_id(census, trim=True)
    census_tag = '.census_only'
plot_tag += census_tag

for t in intron_samples:
    ### sort everything by geneID
    intron_samples[t].sort_by_gene_id()

    ### binarize intron counts into oberserved / not observed
    intron_samples[t].binarize_intron_matrix(options.thresholds[t])

    ### clean up data from samples with no intron expression
    ### remove samples that express 5 or fewer introns
    intron_samples[t].filter_strains_on_mincount(options.min_introns_per_sample)

### subset to introns that are only scarcely seen in normals
if options.max_norm_expression_frac is not None:
    ### look at all samples together
    if not type(options.max_norm_expression_frac) is dict:
        tmp = sp.hstack([intron_samples[x].introns_bin for x in options.normal_set])
        exp_frac = sp.sum(tmp, axis=1).astype('float') / tmp.shape[1]
        del tmp
        k_idx = sp.where(exp_frac <= options.max_norm_expression_frac)[0]
        print >> sys.stdout, 'keep %i of %i introns that are expressed below the given threshold in normals' % (k_idx.shape[0], exp_frac.shape[0])
    ### compute fraction per normal type and take max
    else:
        exp_frac = sp.vstack([sp.mean(intron_samples[x].introns_bin.astype('float'), axis=1) <= options.max_norm_expression_frac[x] for x in options.normal_set])
        k_idx = sp.where(exp_frac.min(axis=0))[0]
        print >> sys.stdout, 'keep %i of %i introns that are expressed below the given threshold in normals' % (k_idx.shape[0], exp_frac.shape[1])

    for t in intron_samples:
        intron_samples[t].reorder_by_row_idx(k_idx)

### filter for introns that have a higher mean expression in normals than in tumors for the
### samples that show the intron
if options.filter_high_normal_expression:
    ### normalize expression counts on remaining samples
    exp_sample.make_unique_by_strain()
    exp_sample.subset_by_strain(sp.hstack([intron_samples[x].strains for x in intron_samples]))
    ### build tumor expression object or normalization
    exp_obj = ExpressionData(exp_sample.cnt.astype('float'), exp_sample.gids, exp_sample.sids)
    exp_obj.libsizenorm(replace_inf=0, replace_nan=0)

    ### replicate gene expression for each intron to appropriatly filter
    (aa, _) = sp.where(intron_samples['tc'].gid == exp_sample.gids[:, sp.newaxis])
    tmp_idx = sp.argsort(_)
    aa = aa[tmp_idx]
    assert(sp.all(intron_samples['tc'].gid == exp_sample.gids[aa]))

    (_, bb) = sp.where(sp.hstack([intron_samples[x].strains for x in sorted(intron_samples)])[:, sp.newaxis] == exp_sample.sids)
    assert(sp.all(sp.hstack([intron_samples[x].strains for x in sorted(intron_samples)]) == exp_sample.sids[bb]))
    tmp_exp = exp_obj.Y[aa, :][:, bb]
    t_mean = sp.sum(tmp_exp[:, exp_sample.is_tumor[bb]] * intron_samples['tc'].introns_bin, axis=1) / sp.c_[sp.sum(intron_samples['tc'].introns_bin, axis=1), sp.ones((intron_samples['tc'].introns_bin.shape[0], ))].max(axis=1)
    n_mean = sp.mean(tmp_exp[:, ~exp_sample.is_tumor[bb]], axis=1)
    tt_idx = sp.where(n_mean > t_mean)[0]

    ### clean-up introns from genes that have lower mean expression in normals
    for t in intron_samples:
        intron_samples[t].reorder_by_row_idx(tt_idx)

### remove introns that are not confirmed in any sample
k_idx = sp.where(sp.sum(sp.r_[[(sp.sum(intron_samples[x].introns_bin, axis=1) > 0) for x in intron_samples]], axis=0) > 0)[0]
for t in intron_samples:
    intron_samples[t].reorder_by_row_idx(k_idx)

### compute overview matrix - average occurrence per cancer type
for t in intron_samples:
    intron_samples[t].add_overview()

print >> sys.stdout, 'Analysing %i introns in %i genes' % (intron_samples['tc'].eid.shape[0], sp.unique(intron_samples['tc'].gid).shape[0])
print >> sys.stdout, 'Confirmed in %s: %i introns' % (','.join(intron_samples.keys()), sp.sum(sp.sum(sp.hstack([intron_samples[x].introns_bin for x in sorted(intron_samples)]), axis=1) > 0))

###
### RANKING
###
### sort introns using the respective strategy

### compute sort index as maximum over cancer types
if options.strategy == 'per_CT':
    tmp = sp.zeros((intron_samples['tc'].introns_bin.shape[0], intron_samples['tc'].types_u.shape[0]))
    for i, ct in enumerate(sorted(intron_samples['tc'].types_u)):
        tmp[:, i] = intron_samples['tc'].ov[:, i] / (1 + intron_samples['tn'].ov[:, i])
    sortval = tmp.max(axis=1)
### compute sort index over all non-tcga samples
elif options.strategy == 'non_tumor':
    tmp = sp.zeros((intron_samples['tc'].introns_bin.shape[0], len(options.normal_set)))
    for i, label in enumerate(options.normal_set):
        tmp[:, i] = sp.mean(intron_samples['tc'].introns_bin, axis=1) / (0.000001 + sp.mean(intron_samples[label].introns_bin, axis=1))
    sortval = tmp.min(axis=1)
### compue sort index as sum of tumor occurrences over sum of normal ocurrences
elif options.strategy == 'sum_non_tumor':
    sortval = (sp.sum(sp.hstack([intron_samples[x].introns_bin for x in options.tumor_set]), axis=1) + 1) \
            / (sp.sum(sp.hstack([intron_samples[x].introns_bin for x in options.normal_set]), axis=1) + 1)
### compute sort index as mean of tumor ocurrences over mean of normal ocurrences 
elif options.strategy == 'mean_non_tumor':
    sortval = sp.mean(sp.hstack([intron_samples[x].introns_bin for x in options.tumor_set]), axis=1) \
           / (sp.mean(sp.hstack([intron_samples[x].introns_bin for x in options.normal_set]), axis=1) + 0.0000001)
### compute sort index as max fisher p-value of all  
elif options.strategy == 'fisher_CT':
    tmp_tu_wi = sp.sum(sp.hstack([intron_samples[x].introns_bin for x in options.tumor_set]), axis=1)   # number of ocurrences in tumor
    tmp_tu_wo = sum([intron_samples[x].introns_bin.shape[1] for x in options.tumor_set]) - tmp_tu_wi    # number of non-occurences in tumor
    tmp_sort = sp.zeros((intron_samples[options.tumor_set[0]].introns_bin.shape[0], 0))
    test_dict = {'BLDR':'BLCA',
                 'BRST':'BRCA',
                 'KIDN':'KIRC',
                 'LUNG':'LUAD',
                 'LUNG':'LUSC'}
    for norm in options.normal_set:
        if norm != 'gt':
            tmp_no_wi = sp.sum(intron_samples[norm].introns_bin, axis=1)        # number of ocurrences in current normal
            tmp_no_wo = intron_samples[norm].introns_bin.shape[1] - tmp_no_wi   # number of non-ocurrences in current normal
            tmp_sort = sp.c_[tmp_sort, sp.array([fisher.pvalue(tmp_no_wi[i], tmp_no_wo[i], tmp_tu_wi[i], tmp_tu_wo[i]).left_tail \
                                                 for i in range(tmp_no_wi.shape[0])])]
    sortval = 1.0 - tmp_sort.max(axis=1)
### 
elif options.strategy == 'fisher_subset':
    tmp_tu_wi = sp.sum(sp.hstack([intron_samples[x].introns_bin for x in options.tumor_set]), axis=1)   # number of ocurrences in tumor
    tmp_tu_wo = sum([intron_samples[x].introns_bin.shape[1] for x in options.tumor_set]) - tmp_tu_wi    # number of non-occurences in tumor
    tmp_sort = sp.zeros((intron_samples[options.tumor_set[0]].introns_bin.shape[0], 0))
    for norm in options.normal_set:
        idx = sp.sort(sp.hstack([sp.where(intron_samples[norm].types == x)[0] for x in fisher_subset[norm]])) 
        tmp_no_wi = sp.sum(intron_samples[norm].introns_bin[:, idx], axis=1)    # number of ocurrences in current normal
        tmp_no_wo = idx.shape[0] - tmp_no_wi                                    # number of non-ocurrences in current normal
        tmp_sort = sp.c_[tmp_sort, sp.array([fisher.pvalue(tmp_no_wi[i], tmp_no_wo[i], tmp_tu_wi[i], tmp_tu_wo[i]).left_tail \
                                             for i in range(tmp_no_wi.shape[0])])]
    sortval = 1.0 - tmp_sort.max(axis=1)
###
elif options.strategy == 'fisher':
    tmp_tu_wi = sp.sum(sp.hstack([intron_samples[x].introns_bin for x in options.tumor_set]), axis=1)   # number of ocurrences in tumor
    tmp_tu_wo = sum([intron_samples[x].introns_bin.shape[1] for x in options.tumor_set]) - tmp_tu_wi    # number of non-occurences in tumor
    tmp_sort = sp.zeros((introns[options.tumor_set[0]].shape[0], 0))
    tmp_enrich = sp.zeros((introns[options.tumor_set[0]].shape[0], 0))
    for norm in options.normal_set:
        tmp_no_wi = sp.sum(intron_samples[norm].introns_bin, axis=1)            # number of ocurrences in current normal
        tmp_no_wo = intron_samples[norm].introns_bin.shape[0] - tmp_no_wi       # number of non-ocurrences in current normal
        tmp_sort = sp.c_[tmp_sort, sp.array([fisher.pvalue(tmp_no_wi[i], tmp_no_wo[i], tmp_tu_wi[i], tmp_tu_wo[i]).left_tail \
                                             for i in range(tmp_no_wi.shape[0])])]
        tmp_enrich = sp.c_[tmp_enrich, ((tmp_tu_wi / (tmp_tu_wi + tmp_tu_wo.astype('float'))) + 1) /\
                                       ((tmp_no_wi / (tmp_no_wi + tmp_no_wo.astype('float'))) + 1)]
    sortval = 1.0 - tmp_sort.max(axis=1)


s1_idx = sp.argsort(sortval)[::-1]
sortval = sortval[s1_idx]


### 
### PLOTTING
###

plot_list = ['tc', 'tn', 'gt']
tick_list = [intron_samples[x].idx_u for x in plot_list]
label_list = []
for t in plot_list:
    label_list.append(['%s (N=%i)' % (x, sp.sum(intron_samples[t].types == x)) for x in intron_samples[t].types_u])
width_list = [options.width, options.width, options.width, 1]

### plot binary matrix of all top event ocurrences with fixed width per sample
axes = []
fig = plt.figure(figsize=(sum(width_list), 8), dpi=160)
gs = gridspec.GridSpec(1, len(plot_list) + 1, width_ratios=width_list)
for i,p in enumerate(plot_list):
    axes.append(fig.add_subplot(gs[0, i]))
    axes[-1].matshow(intron_samples[p].introns_bin[s1_idx[:options.top], :], aspect='auto', cmap=cm.Blues)

for i, ax in enumerate(axes):
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if i == 0:
        axes[i].set_yticks(sp.arange(options.top))
        axes[i].set_yticklabels(intron_samples['tc'].gid[s1_idx[:options.top]])
    axes[i].set_xticks(tick_list[i])
    axes[i].set_xticklabels(label_list[i], rotation=90)
plt.tight_layout()
if len(event_types) > 1:
    plt.savefig('%s/compare_splice%s.all_events.pdf' % (options.plotdir, plot_tag), format='pdf')
    plt.savefig('%s/compare_splice%s.all_events.png' % (options.plotdir, plot_tag), format='png')
else:
    plt.savefig('%s/compare_splice%s.%s.pdf' % (options.plotdir, plot_tag, event_types[0]), format='pdf')
    plt.savefig('%s/compare_splice%s.%s.png' % (options.plotdir, plot_tag, event_types[0]), format='png')

### plot binary matrix of all top event ocurrences with adapted scaling of widths
if False:
    axes = []
    fig = plt.figure(figsize=(sum(width_list) * 2, 80), dpi=160)
    width_ratios = [intron_samples[x].shape[1] for x in plotlist]
    width_ratios.append(1)
    gs = gridspec.GridSpec(1, len(plot_list) + 1, width_ratios=width_ratios)
    for i,p in enumerate(plot_list):
        axes.append(fig.add_subplot(gs[0, i]))
        axes[-1].matshow(intron_samples[p].introns_bin[s1_idx, :], aspect='auto', cmap=cm.Blues)

    for i, ax in enumerate(axes):
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        axes[i].set_xticks(tick_list[i])
        axes[i].set_xticklabels(label_list[i], rotation=90)
    plt.tight_layout()
    if len(event_types) > 1:
        plt.savefig('%s/compare_splice%s.all_events.full.pdf' % (options.plotdir, plot_tag), format='pdf')
        plt.savefig('%s/compare_splice%s.all_events.full.png' % (options.plotdir, plot_tag), format='png')
    else:
        plt.savefig('%s/compare_splice%s.%s.full.pdf' % (options.plotdir, plot_tag, event_types[0]), format='pdf')
        plt.savefig('%s/compare_splice%s.%s.full.png' % (options.plotdir, plot_tag, event_types[0]), format='png')

### plot average of all sample groups on fixed widths
tick_list = [sp.arange(len(intron_samples['tc'].types_u)), sp.arange(len(intron_samples['tc'].types_u)), sp.arange(len(intron_samples['gt'].types_u)), [0]]
delta = 0
axes = []
ax_exclude = []
fig = plt.figure(figsize=(sum(width_list), 8), dpi=160)
gs = gridspec.GridSpec(1, len(plot_list) + 1, width_ratios=width_list)
for i,p in enumerate(plot_list):
    axes.append(fig.add_subplot(gs[0, i]))
    if len(intron_samples[p].ov.shape) < 2:
        cax = axes[-1].matshow(intron_samples[p].ov[s1_idx[:options.top]][:, sp.newaxis], aspect='auto', vmin=0, vmax=1, cmap=cm.Blues)
    else:
        cax = axes[-1].matshow(intron_samples[p].ov[s1_idx[:options.top], :], aspect='auto', vmin=0, vmax=1, cmap=cm.Blues)

tick_list = [sp.arange(len(intron_samples['tc'].types_u)), sp.arange(len(intron_samples['tc'].types_u)), sp.arange(len(intron_samples['gt'].types_u)), [0]]
### set axes labels
for i, ax in enumerate(axes):
    if i in ax_exclude:
        continue
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if i == 0:
        axes[i].set_yticks(sp.arange(options.top))
        axes[i].set_yticklabels(intron_samples[options.tumor_set[0]].gid[s1_idx[:options.top]])
    axes[i].set_xticks(tick_list[i])
    axes[i].set_xticklabels(label_list[i], rotation=90)

plt.tight_layout()
if len(event_types) > 1:
    plt.savefig('%s/compare_splice%s.all_events.overview.pdf' % (options.plotdir, plot_tag), format='pdf')
    plt.savefig('%s/compare_splice%s.all_events.overview.png' % (options.plotdir, plot_tag), format='png')
else:
    plt.savefig('%s/compare_splice%s.%s.overview.pdf' % (options.plotdir, plot_tag, event_types[0]), format='pdf')
    plt.savefig('%s/compare_splice%s.%s.overview.png' % (options.plotdir, plot_tag, event_types[0]), format='png')

ts = options.tumor_set[0]
### report all gene labels sorted by sortval to a file for functional analysis
### if a label occurs multiple times, the highest ranked one is kept
if False:
    if len(event_types) > 1:
        outfile = gzip.open('%s/compare_splice%s.all_events.ranked_genes.txt.gz' % (options.plotdir, plot_tag), 'w')
    else:
        outfile = gzip.open('%s/compare_splice%s.%s.ranked_genes.txt.gz' % (options.plotdir, plot_tag, event_types[0]), 'w')

    _, g_idx = sp.unique(intron_samples[ts].gid[s1_idx], return_index=True)
    print >> outfile, 'gene_id\tgene_name\tevent_id\t' + ':'.join(intron_samples[ts].strains)

    for x in sorted(g_idx):
        print >> outfile, intron_samples[ts].gid[s1_idx][x] + '\t' + \
                          names.get_ID(intron_samples[ts].gid[s1_idx][x].split('.')[0], lookup) + '\t' + \
                          intron_samples[ts].eid[s1_idx][x] + '\t' + ':'.join(intron_samples[ts].introns[s1_idx[x], :].astype('str'))
    outfile.close()

### collapse introns that belong to same gene
labels = sp.array([x.split('.')[0] for x in intron_samples[ts].gid[s1_idx]])
### sort labels
idx = sp.argsort(labels)
labels = labels[idx]
sortval_resorted = sortval[idx]
for t in intron_samples:
    if len(intron_samples[t].ov.shape) > 1:
        intron_samples[t].ov = intron_samples[t].ov[s1_idx[idx], :]
    else:
        intron_samples[t].ov = intron_samples[t].ov[s1_idx[idx]]

### make unique over genes
labels_u, f_idx = sp.unique(labels, return_index=True)
l_idx = sp.r_[f_idx[1:,], labels.shape[0]]
unique_count = labels_u.shape[0]
for t in intron_samples:
    intron_samples[t].add_overview_per_gene(unique_count)
sortval_ = sp.zeros((unique_count, ))
for m in xrange(unique_count):
    max_idx = sortval_resorted[f_idx[m]:l_idx[m]].argmax()
    for t in intron_samples:
        intron_samples[t].update_overview_per_gene(m, f_idx[m] + max_idx)
    sortval_[m] = sortval_resorted[f_idx[m] + max_idx]

### re-sort by sortval_
idx = sp.argsort(sortval_)[::-1]
labels_u = labels_u[idx]
for t in intron_samples:
    intron_samples[t].ov_per_gene = intron_samples[t].ov_per_gene[idx]
sortval_ = sortval_[idx]

### re-set labels according to sort order
labels = sp.array([x.split('.')[0] for x in intron_samples[ts].gid[s1_idx]])

### plot overview
vm = intron_samples['tc'].ov_per_gene.max()
axes = []
fig = plt.figure(figsize=(sum(width_list), 8), dpi=160)
ax_exclude = []

for i,p in enumerate(plot_list):
    axes.append(fig.add_subplot(gs[0, i]))
    if len(intron_samples[p].ov_per_gene.shape) < 2:
        cax = axes[-1].matshow(sp.sqrt(intron_samples[p].ov_per_gene[:options.top, sp.newaxis]), aspect='auto', vmin=0, vmax=sp.sqrt(vm), cmap=cm.Blues)
    else:
        cax = axes[-1].matshow(sp.sqrt(intron_samples[p].ov_per_gene[:options.top, :]), aspect='auto', vmin=0, vmax=sp.sqrt(vm), cmap=cm.Blues)
axes.append(fig.add_subplot(gs[0, len(plot_list)]))
cb = matplotlib.colorbar.ColorbarBase(axes[-1], norm=cax.norm, cmap=cax.cmap, orientation='vertical')
cb.set_label('fraction of samples', fontsize=14) 
axes[-1].set_yticklabels(['%.1f%%' % (x * 100) for x in sp.power(sp.linspace(0.0, sp.sqrt(vm), len(axes[-1].get_yticks())), 2)])
axes[-1].yaxis.tick_right()
ax_exclude.append(len(axes) - 1)

### set axes labels
for i, ax in enumerate(axes):
    if i in ax_exclude:
        continue
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if i == 0:
        axes[i].set_yticks(sp.arange(options.top))
        trans_labels = sp.array([names.get_ID(x, lookup) for x in labels_u[:options.top]])
        trans_labels = sp.array([trans_labels[j] + ' *' if j in enrichment_idxs['census'] else trans_labels[j] for j in range(trans_labels.shape[0])])
        axes[i].set_yticklabels(trans_labels)
    axes[i].set_xticks(tick_list[i])
    axes[i].set_xticklabels(label_list[i], rotation=90)

plt.tight_layout()

if len(event_types) > 1:
    plt.savefig('%s/compare_splice%s.all_events.overview.merged.pdf' % (options.plotdir, plot_tag), format='pdf')
    plt.savefig('%s/compare_splice%s.all_events.overview.merged.png' % (options.plotdir, plot_tag), format='png')
else:
    plt.savefig('%s/compare_splice%s.%s.overview.merged.pdf' % (options.plotdir, plot_tag, event_types[0]), format='pdf')
    plt.savefig('%s/compare_splice%s.%s.overview.merged.png' % (options.plotdir, plot_tag, event_types[0]), format='png')

if options.plot_intron_distribution:

    for j in range(options.top):

        fig = plt.figure(figsize=(40, 15), dpi=80)
        gs = gridspec.GridSpec(1, len(intron_samples))

        for tt,t in enumerate(intron_samples):
            ### plot dist
            ax = fig.add_subplot(gs[0, tt])
            dist = []
            for ct in intron_samples[t].types_u:
                cc_idx = sp.where(intron_samples[t].types == ct)[0]
                ### get median quotients of current intron
                ### get all introns of that gene
                g_idx = sp.where(intron_samples[t].gid == intron_samples[t].gid[s1_idx[j]])[0]
                denom = [sp.mean(intron_samples[t].introns[g_idx, l]) for l in cc_idx] 
                dist.append(sp.array([intron_samples[t].introns[s1_idx[j], l] / denom[k] if denom[k] > 0 else 0 for k,l in enumerate(cc_idx)]))
            violin_plot(ax, dist, range(intron_samples[t].types_u.shape[0]), bp=True)
            ax.set_xticklabels(intron_samples[t].types_u, rotation=45)
            ax.set_title('%s: %s - %s' % (t, labels[j], names.get_ID(labels[j])))

        plt.tight_layout()
        plt.savefig('%s/compare_splice%s.top_%i_intron_distributions.pdf' % (options.plotdir, plot_tag, j + 1), format='pdf')   
        plt.savefig('%s/compare_splice%s.top_%i_intron_distributions.png' % (options.plotdir, plot_tag, j + 1), format='png')   
        plt.close(fig)

