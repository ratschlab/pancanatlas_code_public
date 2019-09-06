import sys
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import numpy.random as npr
npr.seed(23)

import scipy as sp
import scipy.stats as spst
import pdb
import h5py
import os
import re
import gzip
import cPickle

sys.path.append('/cluster/home/akahles/git/projects/2013/PanCancerTCGA/rerun2017')
import utils.names as names
from utils.paths import *
import utils.utils as utils
sys.path.append('/cluster/home/akahles/git/tools/python/viz') 
from distribution import violin_plot

### basic settings
PLOTDIR = os.path.join(paths.plotdir, 'as_complexity')
if not os.path.exists(PLOTDIR):
    os.makedirs(PLOTDIR)
CONF = 2

if not os.path.exists(paths.basedir_tss_fix):
    os.makedirs(paths.basedir_tss_fix)

fl_tag = ''
filt_lib = True
if filt_lib:
    fl_tag = '.filtLib'

conf_count = 3
conf_count_glob = 20
conf_tag = '.globsum%i%s.conf%i' % (conf_count_glob, fl_tag, conf_count)
gtex_thresh = 0.01
counts_tcga = []
counts_gtex = []
donor_support_tcga = []
log = False
use_whitelist = True

log_tag = ''
if log:
    log_tag = '.log10'

wl_tag = ''
libsize = paths.libsize
if use_whitelist:
    wl_tag = '.whitelisted'
    libsize = paths.libsize_whitelist

pickle_tag = wl_tag + '.G%s' % str(gtex_thresh) + conf_tag
plot_tag = log_tag + pickle_tag

# get lookup table of gene names
lookup = names.get_lookup_complete() 

### get TCGA type dictionary
(ct_dict, is_tumor_dict) = utils.get_ct_dict_metatable(paths.metadata, style='pancan_rerun18')

### create GTex type dictionary
sample_dict = {'ADIP':'../annotation/sample_lists_gtex/Adipose.txt',
               'BLDR':'../annotation/sample_lists_gtex/Bladder.txt',
               'BRAI':'../annotation/sample_lists_gtex/Brain.txt',
               'BRST':'../annotation/sample_lists_gtex/Breast.txt',
               'CERV':'../annotation/sample_lists_gtex/Cervix.txt',
               'COLN':'../annotation/sample_lists_gtex/Colon.txt',
               'ESPH':'../annotation/sample_lists_gtex/Esophagus.txt',
               'KIDN':'../annotation/sample_lists_gtex/Kidney.txt',
               'LIVR':'../annotation/sample_lists_gtex/Liver.txt',
               'LUNG':'../annotation/sample_lists_gtex/Lung.txt',
               'OVRY':'../annotation/sample_lists_gtex/Ovary.txt',
               'PANC':'../annotation/sample_lists_gtex/Pancreas.txt',
               'PRST':'../annotation/sample_lists_gtex/Prostate.txt',
               'SKIN':'../annotation/sample_lists_gtex/Skin.txt',
               'STOM':'../annotation/sample_lists_gtex/Stomach.txt',
               'THYR':'../annotation/sample_lists_gtex/Thyroid.txt',
               'TSTS':'../annotation/sample_lists_gtex/Testis.txt',
               'UTER':'../annotation/sample_lists_gtex/Uterus.txt',
               'CELL':'../annotation/sample_lists_gtex/Cells.txt'}
gt_dict = utils.get_gt_dict(sample_dict)

### how shall we map the samples onto each other?
gt_tcga_dict = {'BLDR':['BLCA'],
                'BRAI':['LGG', 'GBM'],
                'BRST':['BRCA'],
                'CERV':['CESC'],
                'COLN':['COAD', 'READ'],
                'ESPH':['ESCA'],
                'KIDN':['KICH', 'KIRC', 'KIRP', 'ACC', 'PCPG'],
                'LIVR':['LIHC'],
                'LUNG':['LUAD', 'LUSC'],
                'OVRY':['OV'], 
                'PANC':['PAAD'],
                'PRST':['PRAD'],
                'SKIN':['SKCM'],
                'STOM':['STAD'],
                'THYR':['THCA'],
                'TSTS':['TGCT'],
                'UTER':['UCEC', 'UCS']}

### get list of protein coding genes
cod_genes = sp.loadtxt(paths.coding_genes, dtype='str', delimiter='\t')
cod_genes = sp.array(cod_genes[:, 4], dtype='str')

### get list of unique cancer/tissue types used
ctypes_u = sp.unique(sp.array(ct_dict.values(), dtype='str'))
gtypes_u = sp.unique(sp.array(gt_dict.values(), dtype='str'))

### preprocessing to harmonize genes used in the analysis
print 'loading data for TCGA'
IN_TC = h5py.File(os.path.join(paths.basedir_as, 'spladder', 'genes_graph_conf%i.merge_graphs.validated.count%s.hdf5' % (CONF, wl_tag)), 'r')
gids_tcga = IN_TC['gene_ids_edges'][:, 0]
gnames_tcga = IN_TC['gene_names'][:]
strains_tcga = IN_TC['strains'][:]
strains_tcga_short = sp.array(['-'.join(x.split('-')[:3]) for x in strains_tcga])
ctypes = sp.array([ct_dict[x.split('.')[1]] for x in strains_tcga], dtype='str')
tcga_is_tumor = sp.array([is_tumor_dict[x.split('.')[1]] for x in strains_tcga], dtype='bool')

if filt_lib:
    ### get libsize information
    data_lib = sp.loadtxt(libsize, delimiter='\t', dtype='str')
    data_lib = data_lib[1:, :]
    data_lib[:, 0] = sp.array([x.split('.')[0] for x in data_lib[:, 0]])
    a,b = sp.where(sp.array([x.split('.')[0] for x in strains_tcga])[:, sp.newaxis] == data_lib[:, 0])
    assert sp.all(sp.array([x.split('.')[0] for x in strains_tcga]) == data_lib[b, 0])
    data_lib = data_lib[b, :]

    lkidx = sp.where(data_lib[:, 1].astype('float') > 2500)[0]
    strains_tcga = strains_tcga[lkidx]
    strains_tcga_short = strains_tcga_short[lkidx]
    ctypes = ctypes[lkidx]
    tcga_is_tumor = tcga_is_tumor[lkidx]
else:
    lkidx = sp.arange(strains_tcga.shape[0])

print 'loading data for GTEx'
IN_GT = h5py.File(os.path.join(paths.basedir_as_gtex, 'spladder', 'genes_graph_conf%i.merge_graphs.validated.count.hdf5' % CONF), 'r')
gids_gtex =  IN_GT['gene_ids_edges'][:, 0]
gnames_gtex = IN_GT['gene_names'][:]
strains_gtex = IN_GT['strains'][:]
gtypes = sp.array([gt_dict[x.split('.')[0]] if x.split('.')[0] in gt_dict else 'NA' for x in strains_gtex], dtype='str')

gid_names_tcga = sp.array([gnames_tcga[i] for i in gids_tcga], dtype='str')
gid_names_gtex = sp.array([gnames_gtex[i] for i in gids_gtex], dtype='str')
gid_names_common = sp.intersect1d(gid_names_tcga, gid_names_gtex)
kidx_tcga = sp.where(sp.in1d(gid_names_tcga, gid_names_common))[0]
kidx_gtex = sp.where(sp.in1d(gid_names_gtex, gid_names_common))[0]
gids_tcga = gids_tcga[kidx_tcga]
gids_gtex = gids_gtex[kidx_gtex]

if not os.path.exists(os.path.join(paths.basedir_tss, 'tss_size_factors%s%s.cpickle' % (wl_tag, fl_tag))):

    ### compute total edge count for GTEx samples
    print 'Computing total edge count for GTEx samples'
    ### get gene intervals
    s_idx = sp.argsort(gids_gtex, kind='mergesort')
    _, f_idx = sp.unique(gids_gtex[s_idx], return_index=True)
    l_idx = sp.r_[f_idx[1:], gids_gtex.shape[0]]
    ### get counts
    genecounts_gtex = sp.zeros((f_idx.shape[0], IN_GT['edges'].shape[1]), dtype='int')
    for i in xrange(f_idx.shape[0]):
        if (i + 1) % 20 == 0:
            sys.stdout.write('.')
            if (i + 1) % 1000 == 0:
                sys.stdout.write('%i/%i\n' % (i + 1, f_idx.shape[0]))
            sys.stdout.flush()
        genecounts_gtex[i, :] = sp.sum(IN_GT['edges'][kidx_gtex[s_idx[f_idx[i]:l_idx[i]]], :], axis=0)
    print 'computing size factors for normalization'
    sf_gtex = spst.scoreatpercentile(genecounts_gtex, 75, axis=0)
    sf_gtex = sp.median(sf_gtex) / sf_gtex
    sf_gtex[sp.isnan(sf_gtex) | sp.isinf(sf_gtex)] = 1
    del genecounts_gtex

    ### compute total edge count for TCGA samples
    print 'Computing total edge count for TCGA samples'
    ### get gene intervals
    s_idx = sp.argsort(gids_tcga, kind='mergesort')
    _, f_idx = sp.unique(gids_tcga[s_idx], return_index=True)
    l_idx = sp.r_[f_idx[1:], gids_tcga.shape[0]]
    ### get counts
    genecounts_tcga = sp.zeros((f_idx.shape[0], lkidx.shape[0]), dtype='int')
    for i in xrange(f_idx.shape[0]):
        if (i + 1) % 20 == 0:
            sys.stdout.write('.')
            if (i + 1) % 1000 == 0:
                sys.stdout.write('%i/%i\n' % (i + 1, f_idx.shape[0]))
            sys.stdout.flush()
        genecounts_tcga[i, :] = sp.sum(IN_TC['edges'][kidx_tcga[s_idx[f_idx[i]:l_idx[i]]], :][:, lkidx], axis=0)
    print 'computing size factors for normalization'
    sf_tcga = spst.scoreatpercentile(genecounts_tcga, 75, axis=0)
    sf_tcga = sp.median(sf_tcga) / sf_tcga
    sf_tcga[sp.isnan(sf_tcga) | sp.isinf(sf_tcga)] = 1
    del genecounts_tcga

    ### store results in pickle
    cPickle.dump((sf_gtex, sf_tcga), open(os.path.join(paths.basedir_tss, 'tss_size_factors%s%s.cpickle' % (wl_tag, fl_tag)), 'w'), -1)
else:
    print 'loading size factors from pickle'
    (sf_gtex, sf_tcga) = cPickle.load(open(os.path.join(paths.basedir_tss, 'tss_size_factors%s%s.cpickle' % (wl_tag, fl_tag)), 'r'))

compl_pickle = os.path.join(paths.basedir_tss, 'tss_complexity_counts%s.cpickle' % pickle_tag)
compl_pickle_junc = os.path.join(paths.basedir_tss, 'tss_complexity_counts%s_junc.hdf5' % pickle_tag)
if not os.path.exists(compl_pickle):

    ### remove edges that are confirmed with at least two reads in more than 5% of the GTEx sample
    for i in range(0, IN_GT['edges'].shape[0], 1000):
        sys.stdout.write('.')
        if i > 0 and i % 10000 == 0:
            sys.stdout.write('%i/%i\n' % (i, IN_GT['edges'].shape[0]))
        sys.stdout.flush()
        tmp = sp.mean((IN_GT['edges'][i:i+1000, :] * sf_gtex) >= 2, axis=1)
        if i == 0:
            k_idx = sp.where(tmp < gtex_thresh)[0]
        else:
            k_idx = sp.r_[k_idx, sp.where(tmp < gtex_thresh)[0] + i]
    print 'Removed %i of %i edges that were common in GTEx' % (IN_GT['edges'].shape[0] - k_idx.shape[0], IN_GT['edges'].shape[0])
    print 'Retaining %i edges' % k_idx.shape[0]

    ### remove edges that are not confirmed with at least conf_count_glob many reads in all TCGA samples
    for i in range(0, IN_TC['edges'].shape[0], 1000):
        sys.stdout.write('.')
        if i > 0 and i % 10000 == 0:
            sys.stdout.write('%i/%i\n' % (i, IN_GT['edges'].shape[0]))
        sys.stdout.flush()
        tmp = (IN_TC['edges'][i:i+1000, :][:, lkidx] * sf_tcga).sum(axis=1)
        if i == 0:
            k_idx_ = sp.where(tmp > conf_count_glob)[0]
        else:
            k_idx_ = sp.r_[k_idx_, sp.where(tmp > conf_count_glob)[0] + i]
    kk_idx = sp.in1d(k_idx, k_idx_)
    print 'Removed %i of %i edges that were not supported with at least %i reads in the TCGA cohort' % (k_idx.shape[0] - sp.sum(kk_idx), k_idx.shape[0], conf_count_glob)
    print 'Retaining %i edges' % sp.sum(kk_idx)
    k_idx = k_idx[kk_idx]


    ### remove edges that are annotated
    INJT = h5py.File(os.path.join(paths.basedir_as, 'spladder', 'genes_graph_conf%i.merge_graphs.validated.junctions%s.hdf5' % (CONF, wl_tag)), 'r')
    post = sp.c_[INJT['chrms'][:], INJT['strand'][:], INJT['pos'][:].astype('str')]
    post = sp.array(['.'.join(x) for x in post])
    INJT.close()
    INJA = h5py.File(paths.anno_junctions, 'r')
    posj = sp.c_[INJA['chrms'][:], INJA['strand'][:], INJA['pos'][:].astype('str')]
    posj = sp.array(['.'.join(x) for x in posj])
    INJA.close()

    kk_idx = sp.where(~sp.in1d(post[k_idx], posj))[0]
    print 'Removed %i of %i edges that were found in the annotation' % (k_idx.shape[0] - kk_idx.shape[0], k_idx.shape[0])
    print 'Retaining %i edges' % kk_idx.shape[0]
    k_idx = k_idx[kk_idx]

    ### make sure that we work in the same coordinate system
    assert sp.all(gids_gtex == gids_tcga)

    ### process GTEx samples
    print '\ncollecting counts for GTEx'
    ### get gene intervals
    s_idx = sp.argsort(gids_gtex[k_idx], kind='mergesort')
    _, f_idx = sp.unique(gids_gtex[k_idx][s_idx], return_index=True)
    l_idx = sp.r_[f_idx[1:], k_idx.shape[0]]

    ### get counts
    for i in xrange(f_idx.shape[0]):
        if (i + 1) % 20 == 0:
            sys.stdout.write('.')
            if (i + 1) % 1000 == 0:
                sys.stdout.write('%i/%i\n' % (i + 1, f_idx.shape[0]))
            sys.stdout.flush()
        counts_gtex.append(sp.sum((IN_GT['edges'][sorted(k_idx[s_idx[f_idx[i]:l_idx[i]]]), :] * sf_gtex) >= conf_count, axis=0))

    ### process TCGA samples
    print 'collecting counts for TCGA'
    ### get gene intervals
    s_idx = sp.argsort(gids_tcga[k_idx], kind='mergesort')
    _, f_idx = sp.unique(gids_tcga[k_idx][s_idx], return_index=True)
    l_idx = sp.r_[f_idx[1:], k_idx.shape[0]]

    junc_used = sp.zeros((k_idx.shape[0], lkidx.shape[0]), dtype='bool')

    ### get counts
    for i in xrange(f_idx.shape[0]):
        if (i + 1) % 20 == 0:
            sys.stdout.write('.')
            if (i + 1) % 1000 == 0:
                sys.stdout.write('%i/%i\n' % (i + 1, f_idx.shape[0]))
            sys.stdout.flush()
        tmp = IN_TC['edges'][sorted(k_idx[s_idx[f_idx[i]:l_idx[i]]]), :][:, lkidx] * sf_tcga
        counts_tcga.append(sp.sum(tmp >= conf_count, axis=0))
        donor_support_tcga.extend(sp.sum(tmp >= conf_count, axis=1))
        junc_used[s_idx[f_idx[i]:l_idx[i]], :] = (tmp >= conf_count)

    IN_GT.close()
    IN_TC.close()

    counts_tcga = sp.array(counts_tcga, dtype='int').T
    counts_gtex = sp.array(counts_gtex, dtype='int').T
    donor_support_tcga = sp.array(donor_support_tcga, dtype='int')

    cPickle.dump((counts_tcga, counts_gtex, donor_support_tcga, k_idx), open(compl_pickle, 'w'), -1)
    JU = h5py.File(compl_pickle_junc, 'w')
    JU.create_dataset(name='junc_used', data=junc_used, compression='gzip')
    JU.close()
else:
    print 'loading counts from pickle'
    (counts_tcga, counts_gtex, donor_support_tcga, k_idx) = cPickle.load(open(compl_pickle, 'r'))
    JU = h5py.File(compl_pickle_junc, 'r')
    junc_used = JU['junc_used'][:]
    JU.close()


### write junctions to file
INJT = h5py.File(os.path.join(paths.basedir_as, 'spladder', 'genes_graph_conf%i.merge_graphs.validated.junctions%s.hdf5' % (CONF, wl_tag)), 'r')
jt_pos = INJT['pos'][:]
jt_chrm = INJT['chrms'][:]
jt_strand = INJT['strand'][:]
INJT.close()
junc_out = gzip.open(os.path.join(paths.basedir_tss, 'tss_complexity_counts%s_neojunctions.tsv.gz' % pickle_tag), 'w')
for i in range(junc_used.shape[1]):
    strain = re.sub(r'.aligned', '', strains_tcga[i])
    iidx = sp.where(junc_used[:, i])[0]
    for ii in iidx:
        print >> junc_out,'\t'.join([strain, gid_names_tcga[k_idx[ii]], '%s:%s:%i-%i' % (jt_chrm[k_idx[ii]], jt_strand[k_idx[ii]], jt_pos[k_idx[ii], 0], jt_pos[k_idx[ii], 1])])
junc_out.close()

means_tcga = sp.mean(counts_tcga, axis=1)
medians_tcga = sp.median(counts_tcga, axis=1)
sums_tcga = sp.sum(counts_tcga, axis=1)

tcga_norm_sums_median = sp.median(sums_tcga[~tcga_is_tumor])
tcga_norm_sums_q25 = spst.scoreatpercentile(sums_tcga[~tcga_is_tumor], 25)
tcga_norm_sums_q75 = spst.scoreatpercentile(sums_tcga[~tcga_is_tumor], 75)
tcga_norm_iqr = (tcga_norm_sums_q75 - tcga_norm_sums_q25)

means_gtex = sp.mean(counts_gtex, axis=1)
medians_gtex = sp.median(counts_gtex, axis=1)
sums_gtex = sp.sum(counts_gtex, axis=1)

### print stats
for c in [1, 2, 10, 100, 1000]:
    print 'number of introns confirmed in at least %i donors: %i' % (c, sp.sum(donor_support_tcga >= c)) 
print 'max support of donors: %i' % donor_support_tcga.max()
print 'median support of donors: %i' % sp.median(donor_support_tcga)

### plot per donor support histogram
fig = plt.figure(figsize=(10, 6), dpi=200)
ax = fig.add_subplot(111)
n, bins, patches = ax.hist(donor_support_tcga, range=(1, 1000), bins=100) 
ax.set_ylabel('Donor support per edge')
ax.set_xlabel('Distribution of support')
#ax.set_xticks([])
#ax.set_title('Splicing complexity per gene')
ax.set_xlim([0, ax.get_xlim()[1]])
axs.set_ticks_outer(ax)
axs.clean_axis(ax)
plt.savefig(os.path.join(PLOTDIR, 'donor_edge_support%s.pdf' % (plot_tag)), format='pdf', bbox_inches='tight')
plt.savefig(os.path.join(PLOTDIR, 'donor_edge_support%s.png' % (plot_tag)), format='png', bbox_inches='tight')
plt.close(fig)


### plot distribution of 

### write outlier samples to disk
outlier_idx = sp.where(means_tcga >= 0.4)[0]
data = sp.vstack([sp.array([x.split('.')[0] for x in strains_tcga[outlier_idx]], dtype='str'), ctypes[outlier_idx], sums_tcga[outlier_idx].astype('str')])
data = sp.hstack([sp.array(['sample_id', 'cancer_type', 'splice_degree'])[:, sp.newaxis], data]).T 
sp.savetxt(open(os.path.join(paths.basedir_tss, 'splicing_complexity%s.outliers.tsv' % pickle_tag), 'w'), data, fmt='%s', delimiter='\t')

### write all samples to disk
data = sp.vstack([sp.array([x.split('.')[0] for x in strains_tcga], dtype='str'), ctypes, sums_tcga.astype('str')])
data = sp.hstack([sp.array(['sample_id', 'cancer_type', 'splice_degree'])[:, sp.newaxis], data]).T 
sp.savetxt(open(os.path.join(paths.basedir_tss, 'splicing_complexity%s.tsv' % pickle_tag), 'w'), data, fmt='%s', delimiter='\t')

### write list of genes sorted by mean splicing burden over all samples
s_idx = sp.argsort(gids_tcga[k_idx], kind='mergesort')
_, f_idx = sp.unique(gids_tcga[k_idx][s_idx], return_index=True)
gene_mean_tcga = sp.mean(counts_tcga, axis=0)
gene_idx_tcga = sp.argsort(gene_mean_tcga)[::-1]
gnames_tmp = gnames_tcga[gids_tcga[k_idx][s_idx][f_idx]][gene_idx_tcga]
data = sp.vstack([gnames_tmp, sp.array([names.get_ID(x.split('.')[0], lookup) for x in gnames_tmp]), gene_mean_tcga[gene_idx_tcga].astype('str')]).T
sp.savetxt(open(os.path.join(paths.basedir_tss, 'splicing_complexity_by_gene%s.tsv' % pickle_tag), 'w'), data, fmt='%s', delimiter='\t') 


#
#
# PLOTTING
#
#

### get TCGA colors
colors = tc.get_color_scheme('study', labels=ctypes_u)
color_dict = dict([(ctypes_u[_], colors[_]) for _ in range(len(colors))])

### plot settings
cmap = plt.get_cmap('spring')
norm = plt.Normalize(0, len(gt_tcga_dict))
gtex_grey = '#929292'


#
# plot distribution of all genes by mean splicing complexity
#
fig = plt.figure(figsize=(14, 3), dpi=200)
ax = fig.add_subplot(111)
ax.plot(sp.arange(gene_idx_tcga.shape[0]), gene_mean_tcga[gene_idx_tcga], 'o', color='b', markeredgecolor='none', alpha=0.7)
ax.set_ylabel('Avg. complexity over samples')
ax.set_xlabel('Genes')
ax.set_xticks([])
ax.set_title('Splicing complexity per gene')
axs.set_ticks_outer(ax)
axs.clean_axis(ax)
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_mean_per_gene%s.pdf' % (plot_tag)), format='pdf', bbox_inches='tight')
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_mean_per_gene%s.png' % (plot_tag)), format='png', bbox_inches='tight')
plt.close(fig)

#
# plot distribution of sums per cancer type - violin plot
#
fig = plt.figure(figsize=(18, 9), dpi=100)
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
ax = fig.add_subplot(gs[1, 0])
### plot sums
labels = []
ymax = 0
lines = []
for j,tt in enumerate(sorted(gt_tcga_dict)):
    # plot TCGA 
    for t in gt_tcga_dict[tt]:
        # tumor
        t_idx = sp.where((ctypes == t) & (tcga_is_tumor))[0]
        if t_idx.size != 0:
            violin_plot(ax, [sums_tcga[t_idx]], [len(labels)], bp=False, fc=color_dict[t])
            labels.append(t + '-T (%i)' % t_idx.shape[0])
        # normal 
        t_idx = sp.where((ctypes == t) & (~tcga_is_tumor))[0]
        if t_idx.size != 0:
            violin_plot(ax, [sums_tcga[t_idx]], [len(labels)], bp=False, fc=color_dict[t])
            labels.append(t + '-N (%i)' % t_idx.shape[0])
        ### plot dashed line to distinguish tumor types
        #ax.plot([len(labels) - 0.5, len(labels) - 0.5], [0, 10], '--', color='0.15')
        #ymax = max(ymax, ax.get_ylim()[1])
        if t_idx.shape[0] > 0:
            ymax = max(ymax, sums_tcga[t_idx].max())
    # plot GTEx
    t_idx = sp.where(gtypes == tt)[0]
    violin_plot(ax, [sums_gtex[t_idx]], [len(labels)], bp=False, fc=gtex_grey)
    labels.append(tt + '-G (%i)' % t_idx.shape[0])

    if j < len(gt_tcga_dict) - 1:
        lines.append(len(labels))
    
### plot dashed line to distinguish groups
for j in lines:
    ax.plot([j - 0.5, j - 0.5], [0, ymax], '--', color='grey')
ax.set_ylim([0, ymax])

ax.set_ylabel('Novel junctions')
#ax.set_xlabel('Samples')
ax.set_xticks(sp.arange(len(labels)))
ax.set_xticklabels(labels, rotation=90)
ax.set_xlim([-1, len(labels)])
axs.set_ticks_outer(ax)
axs.clean_axis(ax)
ax.set_title('Splicing Complexity per cancer type')

if log:
    sums_tcga = sp.log10(sums_tcga)
    sums_gtex = sp.log10(sums_gtex)

### plot distribution of sums per cancer type - Getz plot
ax = fig.add_subplot(gs[0, 0])
### plot sums
labels = []
xticks = []
cumx = 0
buffsize = 200
arrowlen = 0.25
for j,tt in enumerate(sorted(gt_tcga_dict)):
    # plot TCGA 
    for t in gt_tcga_dict[tt]:
        # tumor
        t_idx = sp.where(ctypes == t)[0]
        if t_idx.size != 0:
            s_idx = sp.argsort(sums_tcga[t_idx])
            # tumor
            tt_idx = sp.where(tcga_is_tumor[t_idx][s_idx])[0]
            if tt_idx.shape[0] != 0:
                ax.plot(sp.arange(tt_idx.shape[0]) + cumx, sums_tcga[t_idx[s_idx][tt_idx]],'o',  color=color_dict[t], markeredgecolor='none')
            labels.append(t)
            xticks.append(cumx + int(tt_idx.shape[0] / 2))
            cumx += tt_idx.shape[0] + buffsize
            if True:
                # normal
                tn_idx = sp.where(~tcga_is_tumor[t_idx][s_idx])[0]
                if tn_idx.shape[0] != 0:
                    ax.plot(sp.arange(tn_idx.shape[0]) + cumx, sums_tcga[t_idx[s_idx][tn_idx]],'^',  color=color_dict[t], markeredgecolor='none')
                    labels.append(t + '-N')
                    xticks.append(cumx + int(tn_idx.shape[0] / 2))
                cumx += tn_idx.shape[0] + buffsize
if log:
    ax.set_ylabel('Novel junctions (log10)')
else:
    ax.set_ylabel('Novel junctions')
#ax.set_xlabel('Samples')
ax.set_xticks(xticks)
ax.set_xticklabels(labels, rotation=90)
ax.set_xlim([-1 * buffsize, cumx + buffsize])
ax.set_ylim([0, ymax])
axs.set_ticks_outer(ax)
axs.clean_axis(ax)
ax.set_title('Splicing Complexity per cancer type')

plt.tight_layout()
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_sum%s.pdf' % (plot_tag)), format='pdf')
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_sum%s.png' % (plot_tag)), format='png')
plt.close(fig)


### plot distribution of sums per cancer type - Getz plot - SINGLE T vs N
fig = plt.figure(figsize=(18, 5), dpi=100)
ax = fig.add_subplot(111)
### plot sums
labels = []
xticks = []
cumx = 0
buffsize = 200
arrowlen = 0.25
sort_means = []
for ct in ctypes_u:
    tt_idx = sp.where((ctypes == ct) & tcga_is_tumor)[0]
    sort_means.append(sp.mean(sums_tcga[tt_idx]))
sort_means = sp.array(sort_means)
sidx = sp.argsort(sort_means)

for t in ctypes_u[sidx]:
    # tumor
    t_idx = sp.where(ctypes == t)[0]
    if t_idx.size != 0:
        s_idx = sp.argsort(sums_tcga[t_idx])
        # tumor
        tt_idx = sp.where(tcga_is_tumor[t_idx][s_idx])[0]
        if tt_idx.size != 0:
            ax.plot(sp.arange(tt_idx.shape[0]) + cumx, sums_tcga[t_idx[s_idx][tt_idx]],'o',  color=color_dict[t], markeredgecolor='none')
            labels.append(t + '-T')
            xticks.append(cumx + int(tt_idx.shape[0] / 2))
        cumx += tt_idx.shape[0] + buffsize

        # normal
        tn_idx = sp.where(~tcga_is_tumor[t_idx][s_idx])[0]
        if tn_idx.size != 0:
            ax.plot(sp.arange(tn_idx.shape[0]) + cumx, sums_tcga[t_idx[s_idx][tn_idx]],'^',  color=color_dict[t], markeredgecolor='none')
            labels.append(t + '-N')
            xticks.append(cumx + int(tn_idx.shape[0] / 2))
        cumx += tn_idx.shape[0] + buffsize

if log:
    ax.set_ylabel('Novel junctions (log10)')
else:
    ax.set_ylabel('Novel junctions')
#ax.set_xlabel('Samples'
ax.set_xticks(xticks)
ax.set_xticklabels(labels, rotation=90)
ax.set_xlim([-1 * buffsize, cumx + buffsize])
ax.set_ylim([0, ymax])
ax.grid(b=True, which='major', linestyle='--', linewidth=0.2, color='#222222')
ax.xaxis.grid(False)
axs.set_ticks_outer(ax)
axs.clean_axis(ax)
ax.set_title('Splicing Complexity per cancer type')

plt.tight_layout()
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_sum_single_TN%s.pdf' % (plot_tag)), format='pdf')
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_sum_single_TN%s.png' % (plot_tag)), format='png')
plt.close(fig)





### plot distribution of sums per cancer type - Getz plot - SINGLE
fig = plt.figure(figsize=(18, 5), dpi=100)
ax = fig.add_subplot(111)
### plot sums
labels = []
xticks = []
cumx = 0
buffsize = 200
arrowlen = 0.25
sort_means = []
for ct in ctypes_u:
    tt_idx = sp.where((ctypes == ct) & tcga_is_tumor)[0]
    sort_means.append(sp.mean(sums_tcga[tt_idx]))
sort_means = sp.array(sort_means)
sidx = sp.argsort(sort_means)

for t in ctypes_u[sidx]:
    # tumor
    t_idx = sp.where(ctypes == t)[0]
    if t_idx.size != 0:
        s_idx = sp.argsort(sums_tcga[t_idx])
        # tumor
        tt_idx = sp.where(tcga_is_tumor[t_idx][s_idx])[0]
        if tt_idx.shape[0] != 0:
            ax.plot(sp.arange(tt_idx.shape[0]) + cumx, sums_tcga[t_idx[s_idx][tt_idx]],'o',  color=color_dict[t], markeredgecolor='none')
            tn_idx = sp.where(~tcga_is_tumor[t_idx][s_idx])[0]
            if tn_idx.shape[0] >= 5:
                tn_median = sp.median(sums_tcga[t_idx[s_idx][tn_idx]])
                ax.plot([cumx-25, t_idx.shape[0]+cumx+25], [tn_median, tn_median], ':r', linewidth=2.0)
        labels.append(t)
        xticks.append(cumx + int(t_idx.shape[0] / 2))
        cumx += t_idx.size + buffsize
print 'Average number of junctions over tumor types: %.2f' % sp.mean(sums_tcga[tcga_is_tumor])
print 'Min number of junctions over tumor types: %.2f' % sums_tcga[tcga_is_tumor].min()
print 'Max number of junctions over tumor types: %.2f' % sums_tcga[tcga_is_tumor].max()

if log:
    ax.set_ylabel('Novel junctions (log10)')
else:
    ax.set_ylabel('Novel junctions')
#ax.set_xlabel('Samples'
ax.set_xticks(xticks)
ax.set_xticklabels(labels, rotation=90)
ax.set_xlim([-1 * buffsize, cumx + buffsize])
ax.set_ylim([0, ymax])
ax.grid(b=True, which='major', linestyle='--', linewidth=0.2, color='#222222')
ax.xaxis.grid(False)
axs.set_ticks_outer(ax)
axs.clean_axis(ax)
ax.set_title('Splicing Complexity per cancer type')

plt.tight_layout()
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_sum_single%s.pdf' % (plot_tag)), format='pdf')
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_sum_single%s.png' % (plot_tag)), format='png')
plt.close(fig)






### plot distribution of means per cancer type - Getz plot, sorted by max
buffsize = 200
arrowlen = 0.25
### get max values to determine sort order
sort_vals = []
ctypes_un = sp.unique(ctypes[~tcga_is_tumor])
for ct in ctypes_un:
    #t_idx = sp.where(ctypes == ct)[0]
    #sort_vals.append(means_tcga[t_idx].max())
    t1_idx = sp.where((ctypes == ct) & tcga_is_tumor)[0]
    t2_idx = sp.where((ctypes == ct) & ~tcga_is_tumor)[0]
    sort_vals.append(means_tcga[t1_idx].mean() - means_tcga[t2_idx].mean())
ct_idx = sp.argsort(sort_vals)[::-1]
norm = plt.Normalize(0, len(ctypes_un))
half = len(ctypes_un) / 2

fig = plt.figure(figsize=(12, 6), dpi=100)
gs = gridspec.GridSpec(2, 1)
ax = fig.add_subplot(gs[0, 0])
labels = []
xticks = []
cumx = 0
for j,t in enumerate(ctypes_un[ct_idx]):
    if j > half:
        break
    # plot TCGA 
    if sp.sum((ctypes == t) & ~tcga_is_tumor) < 5:
        continue
    # tumor
    t_idx = sp.where(ctypes == t)[0]
    if t_idx.size != 0:
        s_idx = sp.argsort(means_tcga[t_idx])
        # tumor
        tt_idx = sp.where(tcga_is_tumor[t_idx][s_idx])[0]
        if tt_idx.size != 0:
            ax.plot(tt_idx + cumx, means_tcga[t_idx[s_idx][tt_idx]],'o',  color=color_dict[t], markeredgecolor='none')
        # normal 
        tn_idx = sp.where(~tcga_is_tumor[t_idx][s_idx])[0]
        if tn_idx.size != 0:
            ax.plot(tn_idx + cumx, means_tcga[t_idx[s_idx][tn_idx]],'o',  color=color_dict[t], markeredgecolor='magenta', markeredgewidth=2.0)
        # ticks
        labels.append(t)
        xticks.append(cumx + int(t_idx.shape[0] / 2))
        cumx += t_idx.size + buffsize
if log:
    ax.set_ylabel('Mean number (log10) of junctions per gene')
else:
    ax.set_ylabel('Mean number of junctions per gene')
ax.set_xlabel('Samples')
ax.set_xticks(xticks)
ax.set_xticklabels(labels, rotation=90)
ax.set_xlim([-1 * buffsize, cumx + buffsize])
ax.set_ylim([0, 1])
ax.set_title('Degree of Splicing Abnormality')
axs.set_ticks_outer(ax)
axs.clean_axis(ax)

ax = fig.add_subplot(gs[1, 0])
labels = []
xticks = []
cumx = 0
for j,t in enumerate(ctypes_un[ct_idx]):
    if j <= half:
        continue
    # plot TCGA 
    # tumor
    t_idx = sp.where(ctypes == t)[0]
    if t_idx.size != 0:
        s_idx = sp.argsort(means_tcga[t_idx])
        # tumor
        tt_idx = sp.where(tcga_is_tumor[t_idx][s_idx])[0]
        if tt_idx.size != 0:
            ax.plot(tt_idx + cumx, means_tcga[t_idx[s_idx][tt_idx]],'o',  color=color_dict[t], markeredgecolor='none')
        # normal 
        tn_idx = sp.where(~tcga_is_tumor[t_idx][s_idx])[0]
        if tn_idx.size != 0:
            ax.plot(tn_idx + cumx, means_tcga[t_idx[s_idx][tn_idx]],'o',  color=color_dict[t], markeredgecolor='magenta', markeredgewidth=2.0)
        # ticks
        labels.append(t)
        xticks.append(cumx + int(t_idx.shape[0] / 2))
        cumx += t_idx.size + buffsize

if log:
    ax.set_ylabel('Mean number (log10) of junctions per gene')
else:
    ax.set_ylabel('Mean number of junctions per gene')
ax.set_xlabel('Samples')
ax.set_xticks(xticks)
ax.set_xticklabels(labels, rotation=90)
ax.set_xlim([-1 * buffsize, cumx + buffsize])
ax.set_ylim([0, 1])
axs.set_ticks_outer(ax)
axs.clean_axis(ax)
#ax.set_title('Splicing Complexity per cancer type')

plt.tight_layout()
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_mean%s_max_sorted.pdf' % (plot_tag)), format='pdf')
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_mean%s_max_sorted.png' % (plot_tag)), format='png')
plt.close(fig)


### plot distribution of sums per cancer type - Getz plot, all cancer types
### get max values to determine sort order
sort_vals = []
ctypes_un = sp.unique(ctypes[tcga_is_tumor])
for ct in ctypes_un:
    t_idx = sp.where(ctypes == ct)[0]
    sort_vals.append(sums_tcga[t_idx].max())
ct_idx = sp.argsort(sort_vals)[::-1]
norm = plt.Normalize(0, len(ctypes_un))
half = len(ctypes_un) / 2

fig = plt.figure(figsize=(12, 6), dpi=100)
gs = gridspec.GridSpec(2, 1)
ax = fig.add_subplot(gs[0, 0])
labels = []
xticks = []
cumx = 0
for j,t in enumerate(ctypes_un[ct_idx]):
    if j > half:
        break
    # plot TCGA 
    # tumor
    t_idx = sp.where(ctypes == t)[0]
    if t_idx.size != 0:
        s_idx = sp.argsort(sums_tcga[t_idx])
        # tumor
        tt_idx = sp.where(tcga_is_tumor[t_idx][s_idx])[0]
        if tt_idx.size != 0:
            ax.plot(tt_idx + cumx, sums_tcga[t_idx[s_idx][tt_idx]],'o',  color=color_dict[t], markeredgecolor='none')
        # normal 
        tn_idx = sp.where(~tcga_is_tumor[t_idx][s_idx])[0]
        if tn_idx.size != 0:
            ax.plot(tn_idx + cumx, sums_tcga[t_idx[s_idx][tn_idx]],'o',  color=color_dict[t], markeredgecolor='none') #green', markeredgewidth=2.0)
        # ticks
        labels.append(t)
        xticks.append(cumx + int(t_idx.shape[0] / 2))
        cumx += t_idx.size + buffsize
if log:
    ax.set_ylabel('Splicing Burden (log 10)')
else:
    ax.set_ylabel('Splicing Burden')
#ax.set_xlabel('Samples')
ax.set_xticks(xticks)
ax.set_xticklabels(labels, rotation=90)
ax.set_xlim([-1 * buffsize, cumx + buffsize])
ax.set_yticks(ax.get_yticks()[::2])
ax.grid(b=True, which='major', linestyle='--', linewidth=0.2, color='#222222')
ax.xaxis.grid(False)
ymax = ax.get_ylim()[1]
axs.set_ticks_outer(ax)
axs.clean_axis(ax)

ax = fig.add_subplot(gs[1, 0])
labels = []
xticks = []
cumx = 0
for j,t in enumerate(ctypes_un[ct_idx]):
    if j <= half:
        continue
    # plot TCGA 
    # tumor
    t_idx = sp.where(ctypes == t)[0]
    if t_idx.size != 0:
        s_idx = sp.argsort(sums_tcga[t_idx])
        # tumor
        tt_idx = sp.where(tcga_is_tumor[t_idx][s_idx])[0]
        if tt_idx.size != 0:
            ax.plot(tt_idx + cumx, sums_tcga[t_idx[s_idx][tt_idx]],'o',  color=color_dict[t], markeredgecolor='none')
        # normal 
        tn_idx = sp.where(~tcga_is_tumor[t_idx][s_idx])[0]
        if tn_idx.size != 0:
            ax.plot(tn_idx + cumx, sums_tcga[t_idx[s_idx][tn_idx]],'o',  color=color_dict[t], markeredgecolor='none') #green', markeredgewidth=2.0)
        # ticks
        labels.append(t)
        xticks.append(cumx + int(t_idx.shape[0] / 2))
        cumx += t_idx.size + buffsize

if log:
    ax.set_ylabel('Splicing Burden (log 10)')
else:
    ax.set_ylabel('Splicing Burden')
ax.set_xlabel('Samples')
ax.set_xticks(xticks)
ax.set_xticklabels(labels, rotation=90)
ax.set_xlim([-1 * buffsize, cumx + buffsize])
#ax.set_ylim([0, 1])
ax.set_ylim([0, ymax])
ax.set_yticks(ax.get_yticks()[::2])
ax.grid(b=True, which='major', linestyle='--', linewidth=0.2, color='#222222')
ax.xaxis.grid(False)
ymax = ax.get_ylim()[1]
axs.set_ticks_outer(ax)
axs.clean_axis(ax)
#ax.set_title('Splicing Complexity per cancer type')

plt.tight_layout()
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_sum%s_max_sorted_all.pdf' % (plot_tag)), format='pdf')
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_sum%s_max_sorted_all.png' % (plot_tag)), format='png')
plt.close(fig)

### plot distribution of sums per cancer type - Getz plot, sorted by difference of means - TWO ROWS
buffsize = 200
arrowlen = 0.25
### get max values to determine sort order
sort_vals = []
ctypes_un = sp.unique(ctypes[~tcga_is_tumor])
for ct in ctypes_un:
    t1_idx = sp.where((ctypes == ct) & tcga_is_tumor)[0]
    t2_idx = sp.where((ctypes == ct) & ~tcga_is_tumor)[0]
    sort_vals.append(sums_tcga[t1_idx].mean() - sums_tcga[t2_idx].mean())
ct_idx = sp.argsort(sort_vals)[::-1]
norm = plt.Normalize(0, len(ctypes_un))
half = len(ctypes_un) / 2

fig = plt.figure(figsize=(12, 6), dpi=100)
gs = gridspec.GridSpec(2, 1)
ax = fig.add_subplot(gs[0, 0])
labels = []
xticks = []
cumx = 0
for j,t in enumerate(ctypes_un[ct_idx]):
    if j > half:
        break
    # plot TCGA 
    # tumor
    t_idx = sp.where(ctypes == t)[0]
    if t_idx.size != 0:
        s_idx = sp.argsort(sums_tcga[t_idx])
        # tumor
        tt_idx = sp.where(tcga_is_tumor[t_idx][s_idx])[0]
        if tt_idx.size != 0:
            ax.plot(tt_idx + cumx, sums_tcga[t_idx[s_idx][tt_idx]],'o',  color=color_dict[t], markeredgecolor='none')
        # normal 
        tn_idx = sp.where(~tcga_is_tumor[t_idx][s_idx])[0]
        if tn_idx.size != 0:
            ax.plot(tn_idx + cumx, sums_tcga[t_idx[s_idx][tn_idx]],'o',  color=color_dict[t], markeredgecolor='green', markeredgewidth=2.0)
        # ticks
        labels.append(t)
        xticks.append(cumx + int(t_idx.shape[0] / 2))
        cumx += t_idx.size + buffsize
if log:
    ax.set_ylabel('Splicing Burden (log 10)')
else:
    ax.set_ylabel('Splicing Burden')
#ax.set_xlabel('Samples')
ax.set_xticks(xticks)
ax.set_xticklabels(labels, rotation=90)
ax.set_xlim([-1 * buffsize, cumx + buffsize])
ax.set_yticks(ax.get_yticks()[::2])
ax.grid(b=True, which='major', linestyle='--', linewidth=0.2, color='#222222')
ax.xaxis.grid(False)
ymax = ax.get_ylim()[1]
axs.set_ticks_outer(ax)
axs.clean_axis(ax)

ax = fig.add_subplot(gs[1, 0])
labels = []
xticks = []
cumx = 0
for j,t in enumerate(ctypes_un[ct_idx]):
    if j <= half:
        continue
    # plot TCGA 
    # tumor
    t_idx = sp.where(ctypes == t)[0]
    if t_idx.size != 0:
        s_idx = sp.argsort(sums_tcga[t_idx])
        # tumor
        tt_idx = sp.where(tcga_is_tumor[t_idx][s_idx])[0]
        if tt_idx.size != 0:
            ax.plot(tt_idx + cumx, sums_tcga[t_idx[s_idx][tt_idx]],'o',  color=color_dict[t], markeredgecolor='none')
        # normal 
        tn_idx = sp.where(~tcga_is_tumor[t_idx][s_idx])[0]
        if tn_idx.size != 0:
            ax.plot(tn_idx + cumx, sums_tcga[t_idx[s_idx][tn_idx]],'o',  color=color_dict[t], markeredgecolor='green', markeredgewidth=2.0)
        # ticks
        labels.append(t)
        xticks.append(cumx + int(t_idx.shape[0] / 2))
        cumx += t_idx.size + buffsize

if log:
    ax.set_ylabel('Splicing Burden (log 10)')
else:
    ax.set_ylabel('Splicing Burden')
ax.set_xlabel('Samples')
ax.set_xticks(xticks)
ax.set_xticklabels(labels, rotation=90)
ax.set_xlim([-1 * buffsize, cumx + buffsize])
#ax.set_ylim([0, 1])
ax.set_ylim([0, ymax])
ax.set_yticks(ax.get_yticks()[::2])
ax.grid(b=True, which='major', linestyle='--', linewidth=0.2, color='#222222')
ax.xaxis.grid(False)
ymax = ax.get_ylim()[1]
axs.set_ticks_outer(ax)
axs.clean_axis(ax)


plt.tight_layout()
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_sum%s_meandiff_sorted.pdf' % (plot_tag)), format='pdf')
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_sum%s_meandiff_sorted.png' % (plot_tag)), format='png')
plt.close(fig)



### plot distribution of sums per cancer type - Getz plot, sorted by difference of means - ONE ROW
fig = plt.figure(figsize=(18, 4), dpi=200)
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0, 0])
labels = []
xticks = []
cumx = 0
for j,t in enumerate(ctypes_un[ct_idx]):
    # plot TCGA 
    # tumor
    t_idx = sp.where(ctypes == t)[0]
    if t_idx.size != 0:
        s_idx = sp.argsort(sums_tcga[t_idx])
        # tumor
        tt_idx = sp.where(tcga_is_tumor[t_idx][s_idx])[0]
        if tt_idx.size != 0:
            ax.plot(tt_idx + cumx, sums_tcga[t_idx[s_idx][tt_idx]],'o',  color=color_dict[t], markeredgecolor='none')
        # normal 
        tn_idx = sp.where(~tcga_is_tumor[t_idx][s_idx])[0]
        if tn_idx.size != 0:
            ax.plot(tn_idx + cumx, sums_tcga[t_idx[s_idx][tn_idx]],'o',  color=color_dict[t], markeredgecolor='green', markeredgewidth=2.0)
        # ticks
        labels.append(t)
        xticks.append(cumx + int(t_idx.shape[0] / 2))
        cumx += t_idx.size + buffsize
if log:
    ax.set_ylabel('Splicing Burden (log 10)')
else:
    ax.set_ylabel('Splicing Burden')
#ax.set_xlabel('Samples')
ax.set_xticks(xticks)
ax.set_xticklabels(labels, rotation=90)
ax.set_xlim([-1 * buffsize, cumx + buffsize])
ax.set_yticks(ax.get_yticks()[::2])
ax.grid(b=True, which='major', linestyle='--', linewidth=0.2, color='#222222')
ax.xaxis.grid(False)
ymax = ax.get_ylim()[1]
axs.set_ticks_outer(ax)
axs.clean_axis(ax)
plt.tight_layout()
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_sum%s_meandiff_sorted_singlerow.pdf' % (plot_tag)), format='pdf')
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_sum%s_meandiff_sorted_singlerow.png' % (plot_tag)), format='png')
plt.close(fig)




### plot distribution of sums per cancer type - Getz plot, sorted by difference of means - extreme outliers only - ONE ROW
fig = plt.figure(figsize=(18, 4), dpi=200)
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0, 0])
labels = []
xticks = []
cumx = 0
cnt_out_t = 0
cnt_out_n = 0
IQR_FACT = 2.0
for j,t in enumerate(ctypes_un[ct_idx]):
    # plot TCGA 
    # tumor
    t_idx = sp.where(ctypes == t)[0]
    if t_idx.size != 0:
        s_idx = sp.argsort(sums_tcga[t_idx])
        # tumor - outliers
        tt_idx = sp.where(tcga_is_tumor[t_idx][s_idx])[0]
        if tt_idx.size != 0:
            tt_med = sp.median(sums_tcga[t_idx[s_idx][tt_idx]])
            # non outliers
            ttt_idx = sp.where(sp.absolute(tt_med - sums_tcga[t_idx[s_idx][tt_idx]]) <= IQR_FACT*tcga_norm_iqr)[0]
            if ttt_idx.size != 0:
                ax.plot(tt_idx[ttt_idx] + cumx, sums_tcga[t_idx[s_idx][tt_idx][ttt_idx]],'o',  color='#E0E0E0', markeredgecolor='none', alpha=0.5)
            # outliers
            ttt_idx = sp.where(sp.absolute(tt_med - sums_tcga[t_idx[s_idx][tt_idx]]) > IQR_FACT*tcga_norm_iqr)[0]
            if ttt_idx.size != 0:
                ax.plot(tt_idx[ttt_idx] + cumx, sums_tcga[t_idx[s_idx][tt_idx][ttt_idx]],'o',  color=color_dict[t], markeredgecolor='none')
                cnt_out_t += ttt_idx.shape[0]
        # normal 
        tn_idx = sp.where(~tcga_is_tumor[t_idx][s_idx])[0]
        if tn_idx.size != 0:
            tn_med = sp.median(sums_tcga[t_idx[s_idx][tn_idx]])
            # non outliers
            ttn_idx = sp.where(sp.absolute(tn_med - sums_tcga[t_idx[s_idx][tn_idx]]) <= IQR_FACT*tcga_norm_iqr)[0]
            if ttn_idx.size != 0:
                ax.plot(tn_idx[ttn_idx] + cumx, sums_tcga[t_idx[s_idx][tn_idx][ttn_idx]],'o',  color='#E0E0E0', markeredgecolor='none', alpha=0.5)
            # outliers
            ttn_idx = sp.where(sp.absolute(tn_med - sums_tcga[t_idx[s_idx][tn_idx]]) > IQR_FACT*tcga_norm_iqr)[0]
            if ttn_idx.size != 0:
                ax.plot(tn_idx[ttn_idx] + cumx, sums_tcga[t_idx[s_idx][tn_idx][ttn_idx]],'o',  color=color_dict[t], markeredgecolor='green', markeredgewidth=2.0)
                cnt_out_n += ttn_idx.shape[0]
        # ticks
        labels.append(t)
        xticks.append(cumx + int(t_idx.shape[0] / 2))
        cumx += t_idx.size + buffsize
if log:
    ax.set_ylabel('Splicing Burden (log 10)')
else:
    ax.set_ylabel('Splicing Burden')
#ax.set_xlabel('Samples')
ax.set_xticks(xticks)
ax.set_xticklabels(labels, rotation=90)
ax.set_xlim([-1 * buffsize, cumx + buffsize])
ax.set_yticks(ax.get_yticks()[::2])
ax.grid(b=True, which='major', linestyle='--', linewidth=0.2, color='#222222')
ax.xaxis.grid(False)
ymax = ax.get_ylim()[1]
axs.set_ticks_outer(ax)
axs.clean_axis(ax)
plt.tight_layout()
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_sum%s_meandiff_sorted_outliers_iqr%.1f_only_singlerow.pdf' % (plot_tag, IQR_FACT)), format='pdf')
plt.savefig(os.path.join(PLOTDIR, 'complexity_egde_count_sum%s_meandiff_sorted_outliers_iqr%.1f_only_singlerow.png' % (plot_tag, IQR_FACT)), format='png')
plt.close(fig)
print 'Found %i syndeothryptic outliers in  tumor' % cnt_out_t
print 'Found %i syndeothryptic outliers in  normal' % cnt_out_n

