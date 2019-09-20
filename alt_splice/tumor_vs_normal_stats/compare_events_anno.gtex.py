import sys
import os
import pdb
import cPickle
import h5py
import re

import scipy as sp

sys.path.append('/cluster/home/akahles/git/projects/2013/PanCancerTCGA/rerun2017')
import utils.utils as utils
import utils.samples as samples
from utils.paths import paths

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.append('/cluster/home/akahles/git/tools/python/viz') 
from distribution import violin_plot
from axes import *

sys.path.append('/cluster/home/akahles/git/tools/python/oop_libs')
from expression import ExpressionData

sys.path.append('/cluster/home/akahles/git/projects/2013/PanCancerTCGA/colors')
import TCGA_colors as tc

PLOTDIR = os.path.join(paths.plotdir, 'stats')
CONF = 2
if not os.path.exists(PLOTDIR):
    os.makedirs(PLOTDIR)

metadata_style = 'pancan_rerun18'

event_type_dict = {'exon_skip':'exon skip',
                   'alt_3prime':'alternative 3\' site',
                   'alt_5prime':'alternative 5\' site',
                   'intron_retention':'intron retention'}

tn_only = True
tn_tag = ''
if tn_only:
    tn_tag = '.tn'

if len(sys.argv) > 1 and sys.argv[1] != '--force':
    event_types = sys.argv[1:]
else:
    event_types = ['exon_skip', 'alt_3prime', 'alt_5prime', 'intron_retention']

tissues = [('BLCA', 'BLDR'), ('BRCA', 'BRST'), ('COAD', 'COLN'), ('HNSC', 'SKIN'), ('KICH', 'KIDN'), ('KIRC', 'KIDN'), ('KIRP', 'KIDN'), ('LIHC', 'LIVR'), ('LUAD', 'LUNG'), ('LUSC', 'LUNG'), ('PRAD', 'PRST'), ('READ', 'COLN'), ('STAD', 'STOM'), ('THCA', 'THYR'), ('UCEC', 'UTER')]
tissues_plot = [('BRCA', 'BRST'), ('COAD', 'COLN'), ('HNSC', 'SKIN'), ('KIRC', 'KIDN'), ('LIHC', 'LIVR'), ('LUAD', 'LUNG'), ('LUSC', 'LUNG'), ('PRAD', 'PRST'), ('THCA', 'THYR')]
label_dict = {'gt': {'BLCA':'Bladder (GTEx)',
                     'BRCA':'Breast (GTEx)',
                     'COAD':'Colon (GTEx)',
                     'HNSC':'Skin (GTEx)',
                     'KICH':'Kidney (GTEx)',
                     'KIRC':'Kidney (GTEx)',
                     'KIRP':'Kidney (GTEx)',
                     'LIHC':'Liver (GTEx)',
                     'LUAD':'Lung (GTEx)',
                     'LUSC':'Lung (GTEx)',
                     'PRAD':'Prostate (GTEx)',
                     'READ':'Colon (GTEx)',
                     'STAD':'Stomach (GTEx)',
                     'THCA':'Thyroid (GTEx)',
                     'UCEC':'Uterus (GTEx)'},
              'tc': {'BLCA':'BLCA (T)',
                     'BRCA':'BRCA (T)',
                     'COAD':'COAD (T)',
                     'HNSC':'HNSC (T)',
                     'KICH':'KICH (T)',
                     'KIRC':'KIRC (T)',
                     'KIRP':'KIRP (T)',
                     'LIHC':'LIHC (T)',
                     'LUAD':'LUAD (T)',
                     'PRAD':'PRAD (T)',
                     'READ':'READ (T)',
                     'STAD':'STAD (T)',
                     'THCA':'THCA (T)',
                     'UCEC':'UCEC (T)',
                     'LUSC':'LUSC (T)'},
                     #'ALL':'All (Tumr)'},
              'tn': {'BLCA':'BLCA (N)',
                     'BRCA':'BRCA (N)',
                     'COAD':'COAD (N)',
                     'HNSC':'HNSC (N)',
                     'KICH':'KICH (N)',
                     'KIRC':'KIRC (N)',
                     'KIRP':'KIRP (N)',
                     'LIHC':'LIHC (N)',
                     'LUAD':'LUAD (N)',
                     'PRAD':'PRAD (N)',
                     'READ':'READ (N)',
                     'STAD':'STAD (N)',
                     'THCA':'THCA (N)',
                     'UCEC':'UCEC (N)',
                     'LUSC':'LUSC (N)'}}

### get cancer type dictionary
(ct_dict, is_tumor_dict) = utils.get_ct_dict_metatable(paths.metadata, style=metadata_style)

### create GTex type dictionary
sample_dict = {'BLDR':'../annotation/sample_lists_gtex/Bladder.txt', 
               'BRST':'../annotation/sample_lists_gtex/Breast.txt',
               'KIDN':'../annotation/sample_lists_gtex/Kidney.txt',
               'LIVR':'../annotation/sample_lists_gtex/Liver.txt',
               'LUNG':'../annotation/sample_lists_gtex/Lung.txt',
               'PRST':'../annotation/sample_lists_gtex/Prostate.txt',
               'SKIN':'../annotation/sample_lists_gtex/Skin.txt',
               'STOM':'../annotation/sample_lists_gtex/Stomach.txt',
               'THYR':'../annotation/sample_lists_gtex/Thyroid.txt',
               'UTER':'../annotation/sample_lists_gtex/Uterus.txt',
               'COLN':'../annotation/sample_lists_gtex/Colon.txt',
               'CELL':'../annotation/sample_lists_gtex/Cells.txt'}
gt_dict = utils.get_gt_dict(sample_dict)

nan_t = 10 ### this is the number of non-NaN samples per tissue type for an event to be counted
d_psi_t = [0.0, 0.1, 0.3, 0.5]
mr_t = [5, 20, 50]
sf = None

for event_type in event_types:

    print 'Processing %s' % event_type

    ### events are the same for all three datasets
    anno_file = os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.anno_only.pickle' % (event_type, CONF))
    anno_idx = cPickle.load(open(anno_file, 'r'))
    non_anno_idx = sp.where(~anno_idx)[0]
    anno_idx = sp.where(anno_idx)[0]

    if not os.path.exists(os.path.join(paths.basedir_as, 'comparisons', 'pickle')):
        os.makedirs(os.path.join(paths.basedir_as, 'comparisons', 'pickle'))
    count_pickle = os.path.join(paths.basedir_as, 'comparisons', 'pickle', 'sample_psi_counts.%s.conf_%i.pickle' % (event_type, nan_t))
    ### check whether we can load the already pre-processed version of the counts
    if os.path.exists(count_pickle) and not '--force' in sys.argv:
        print 'loading count data from %s' % count_pickle
        (count, tissues, dsets, tids, is_tumor, t_idx) = cPickle.load(open(count_pickle, 'r'))
    else:
        ### list of hdf5 files to take into account
        hdf5_files = {'tc' : os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, CONF)),
                      'tn' : os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, CONF)),
                      'gt' : os.path.join(paths.basedir_as_gtex, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, CONF))}

        hdf5_handles = dict()
        dsets = ['tc', 'tn', 'gt']
        for t in dsets:
            hdf5_handles[t] = h5py.File(hdf5_files[t], 'r')
        sf = dict()
        for t in dsets:
            sf[t] = None

        tids = dict()
        for t in dsets:
            if t == 'gt':
                tids[t] = sp.array([gt_dict[x.split('.')[0]] if x.split('.')[0] in gt_dict else 'NA' for x in hdf5_handles[t]['strains'][:]])
                sf[t] = sp.ones((tids[t].shape[0],), dtype='float')
            else:
                tids[t] = sp.array([ct_dict[x.split('.')[1]] for x in hdf5_handles[t]['strains'][:]])

                ### load libsize
                libsize = sp.loadtxt(paths.libsize, dtype='str', delimiter='\t')
                libsize = libsize[1:, :]
                libsize[:, 0] = sp.array([x.split('.')[1] for x in libsize[:, 0]])
                strains = sp.array([x.split('.')[1] for x in hdf5_handles[t]['strains'][:]])
                a, b = sp.where(strains[:, sp.newaxis] == libsize[:, 0])
                assert sp.all(strains == libsize[b, 0])
                libsize = libsize[b, :]
                sfm = sp.median(libsize[:, 1].astype('float'))
                sf[t] = sfm / libsize[:, 1].astype('float')

        is_tumor = sp.array([is_tumor_dict[x.split('.')[1]] for x in hdf5_handles['tc']['strains'][:]], dtype='bool')

        count = dict()
        for t in dsets:
            for p in tissues:
                count[(t, p)] = dict()

        t_idx = dict()
        for i, p in enumerate(tissues):
            for j, t in enumerate(dsets):
                if p[0] == 'ALL':
                    if t == 'gt':
                        t_idx[(t, p)] = sp.arange(tids[t].shape[0])
                    elif t == 'tc':
                        t_idx[(t, p)] = sp.where(is_tumor)[0]
                    elif t == 'tn':
                        t_idx[(t, p)] = sp.where(~is_tumor)[0]
                else:
                    if t == 'gt':
                        t_idx[(t, p)] = sp.where(tids[t] == p[1])[0]
                    elif t == 'tc':
                        t_idx[(t, p)] = sp.where((tids[t] == p[0]) & is_tumor)[0]
                    elif t == 'tn':
                        t_idx[(t, p)] = sp.where((tids[t] == p[0]) & ~is_tumor)[0]
 
        for j, t in enumerate(dsets):

            print 'computing delta PSI values for %s' % t

            chunksize = hdf5_handles[t]['psi'].chunks[1] * 2 
            for cc, chunk in enumerate(range(0, hdf5_handles[t]['psi'].shape[1], chunksize)):
                if cc > 0 and cc % 10 == 0:                          
                    sys.stderr.write('.')                             
                    if cc % 100 == 0:
                        sys.stderr.write('%i/%i\n' % (cc, hdf5_handles[t]['psi'].shape[1] / chunksize + 1))
                    sys.stderr.flush()  

                tmp_idx = sp.arange(chunk, min(chunk + chunksize, hdf5_handles[t]['psi'].shape[1]))
                tmp_psi_ = hdf5_handles[t]['psi'][:, tmp_idx]
                tmp_iso1 = hdf5_handles[t]['iso1'][:, tmp_idx] * sf[t][:, sp.newaxis]
                tmp_iso2 = hdf5_handles[t]['iso2'][:, tmp_idx] * sf[t][:, sp.newaxis]
                for mr in mr_t:
                    for i, p in enumerate(tissues):
                        tmp_psi = tmp_psi_[t_idx[(t, p)], :]
                        n_idx = sp.c_[tmp_iso1[t_idx[(t, p)], :].max(axis=0), tmp_iso2[t_idx[(t, p)], :].max(axis=0)].min(axis=1) < mr
                        tmp_psi[:, n_idx] = sp.nan
                        idx_nn = ~sp.isnan(tmp_psi)
                        d_psi = sp.nanmax(tmp_psi, axis=0) - sp.nanmin(tmp_psi, axis=0)
                        d_psi[sp.isnan(d_psi)] = 0
                        for dp in d_psi_t:
                            if cc == 0:
                                count[(t, p)][(mr, dp)] = tmp_idx[(sp.sum(idx_nn, axis=0) >= nan_t) & (d_psi >= dp)]
                            else:
                                count[(t, p)][(mr, dp)] = sp.r_[count[(t, p)][(mr, dp)], tmp_idx[(sp.sum(idx_nn, axis=0) >= nan_t) & (d_psi >= dp)]]

        ### store count as pre-processed pickle         
        cPickle.dump((count, tissues, dsets, tids, is_tumor, t_idx), open(count_pickle, 'w'), -1)

    colors = tc.get_color_scheme('study', labels=[x[0] for x in tissues])
    colors = dict([(x[0], colors[i]) for i, x in enumerate(tissues)])

    if tn_only:
        dsets_plot = ['tc', 'tn']
    else:
        dsets_plot = ['tc', 'tn', 'gt']


    ### PLOTTING 
    for mr in mr_t:

        plot_tag = tn_tag + '.mr%i' % mr 

        ###
        ### bar charts (relative thresholds)
        ###
        fig = plt.figure(figsize=(7, 3*len(d_psi_t)), dpi=80)
        gs = gridspec.GridSpec(len(d_psi_t), 1)
        for k, thresh in enumerate(d_psi_t):
            ax = fig.add_subplot(gs[k, 0])
            labels = []
            xloc = []
            buffsize = 1
            buff = buffsize
            j = 0
            dfacts = []
            for p in tissues:
                if not p in tissues_plot:
                    continue
                ### annotated
                counts = []
                curr_xloc = []
                for i, t in enumerate(dsets_plot):
                    if p[0] == 'ALL':
                        t_idx = sp.arange(tids[t].shape[0])
                    else:
                        if t == 'gt':
                            t_idx = sp.where(tids[t] == p[1])[0]
                        else:
                            t_idx = sp.where(tids[t] == p[0])[0]
                    if count[(t, p)] == 0:
                        counts.append(0)
                    else:
                        counts.append(sp.intersect1d(anno_idx, count[(t, p)][thresh]).shape[0])
                    curr_xloc.append(j*len(dsets_plot) + i + buff)
                    if t == 'tn':
                        labels.append(label_dict[t][p[0]] + '\nN=%i' % (sp.sum(~is_tumor[t_idx])))
                    elif t == 'tc':
                        labels.append(label_dict[t][p[0]] + '\nN=%i' % (sp.sum(is_tumor[t_idx])))
                    else:
                        labels.append(label_dict[t][p[0]] + '\nN=%i' % (t_idx.shape[0]))

                ax.bar(curr_xloc, counts, 0.5, color=colors[p[0]])

                ### not annotated
                counts2 = []
                for i, t in enumerate(dsets_plot):
                    if p[0] == 'ALL':
                        t_idx = sp.arange(tids[t].shape[0])
                    else:
                        if t == 'gt':
                            t_idx = sp.where(tids[t] == p[1])[0]
                        else:
                            t_idx = sp.where(tids[t] == p[0])[0]
                    if count[(t, p)] == 0:
                        counts2.append(0)
                    else:
                        counts2.append(sp.intersect1d(non_anno_idx, count[(t, p)][(mr, thresh)]).shape[0])

                ax.bar(curr_xloc, counts2, 0.5, color=colors[p[0]], alpha=0.5, bottom=counts)
                dfact = float(count[('tc',     p)][(mr, thresh)].shape[0]) / float(count[('tn', p)][(mr, thresh)].shape[0]) * 100
                dfacts.append(dfact)
                print 'for mr: %i dpsi: %.1f ct: %s\tT: %i\tN: %i\tT rel to N: %.2f %%' % (mr, thresh, label_dict['tc'][p[0]], count[('tc', p)][(mr, thresh)].shape[0], count[('tn', p)][(mr, thresh)].shape[0], dfact)
                
                buff += buffsize
                j += 1
                xloc.extend(curr_xloc)
            print 'AVG CHANGE: %.2f %%\n' % sp.mean(sp.array(dfacts, dtype='float'))

            ax.set_title('Confirmed %s events (d_psi >%i%%)' % (event_type_dict[event_type], thresh * 100))
            ax.set_ylabel('Confirmed events')
            ax.set_xlabel('Sample group')
            ax.set_xticks(sp.array(xloc) + 0.25)
            if k == len(d_psi_t) - 1:
                ax.set_xticklabels(labels, rotation=90)
            else:
                ax.set_xticklabels([])
            ax.grid(b=True, which='major', linestyle='--', linewidth=0.2, color='#222222')
            ax.xaxis.grid(False)
            set_ticks_outer(ax)
            clean_axis(ax)

        plt.tight_layout()
        plt.savefig(os.path.join(PLOTDIR, 'dpsi_event_overview.%s%s.pdf' % (event_type, plot_tag)), format='pdf')
        plt.savefig(os.path.join(PLOTDIR, 'dpsi_event_overview.%s%s.png' % (event_type, plot_tag)), format='png')
        plt.close(fig)

