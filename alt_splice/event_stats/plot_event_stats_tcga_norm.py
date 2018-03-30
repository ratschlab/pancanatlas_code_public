import sys
import os
import scipy as sp
import cPickle
import h5py

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.append('/cluster/home/akahles/git/projects/2013/PanCancerTCGA/colors')
import TCGA_colors as tc 

sys.path.append('/cluster/home/akahles/git/tools/python/viz/')
import axes as axs

sys.path.append('/cluster/home/akahles/git/projects/2013/PanCancerTCGA/rerun2017')
from utils.paths import paths

PLOTDIR = os.path.join(paths.plotdir, 'stats')
if not os.path.exists(PLOTDIR):
    os.makedirs(PLOTDIR)
CONF = 2

d_psi_t = [0.0, 0.1, 0.3, 0.5]
mr_t = [5, 20, 50]

event_types = ['exon_skip',
               'intron_retention',
               'alt_3prime',
               'alt_5prime', 
               'mutex_exons']

event_dict = {'exon_skip':'Exon Skip',
              'intron_retention':'Intron Retention',
              'alt_3prime':'Alternative 3\' Site',
              'alt_5prime':'Alternative 5\' Site',
              'mutex_exons':'Mutually Exclusive Exons',
              'mult_exon_skip':'Mult. Exon Skip'}


def main():
    figs = dict()
    gss = dict()
    gss['stats'] = gridspec.GridSpec(3, 2) #, wspace=0.0, hspace=0.0)
    gss['stats_log'] = gridspec.GridSpec(3, 2) #, wspace=0.0, hspace=0.0)

    for e, event_type in enumerate(event_types):

        print >> sys.stderr, 'Handling %s' % event_type
        
        ### load annotation index
        is_anno_gtex = cPickle.load(open(os.path.join(paths.basedir_as_gtex, 'merge_graphs_%s_C%i.anno_only.pickle' % (event_type, CONF)), 'r'))
        is_anno_tcga = cPickle.load(open(os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.anno_only.pickle' % (event_type, CONF)), 'r'))

        ### load confident events
        IN = h5py.File(os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, CONF)), 'r')
        idx_conf_tcga = sp.zeros((is_anno_tcga.shape[0],), dtype='bool')
        idx_conf_tcga[IN['conf_idx'][:]] = 1
        IN.close()

        ### for mutex events load cluster idx
        if event_type == 'mutex_exons':
            mutex_clusters = cPickle.load(open(os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.cluster_idx.pickle' % (event_type, CONF)), 'r'))

        ### load psi filtered events
        idx_psi_tcga = cPickle.load(open(os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.counts.psi_filt_per_ct_normalized.pickle' % (event_type, CONF)), 'r'))[1]

        ### get all cancer types
        ctypes = sp.unique([x[1] for x in idx_psi_tcga.keys()])
        colors = sp.array(tc.get_color_scheme('study', labels=ctypes))

        for dp in d_psi_t:
            for mr in mr_t:
                ### get counts
                if  event_type == 'mutex_exons':
                    counts_anno = []
                    counts_new = []
                    for ct in ctypes:
                        if not (mr, ct, dp) in idx_psi_tcga:
                            counts_anno.append(0)
                            counts_new.append(0)
                            continue
                        tmp_idx = sp.zeros((is_anno_tcga.shape[0],), dtype='bool')
                        tmp_idx[idx_psi_tcga[(mr, ct, dp)]] = 1
                        counts_anno.append(sp.unique(mutex_clusters[idx_conf_tcga & is_anno_tcga & tmp_idx]).shape[0])
                        counts_new.append(sp.unique(mutex_clusters[idx_conf_tcga & ~is_anno_tcga & tmp_idx]).shape[0])
                    counts_anno = sp.array(counts_anno)
                    counts_new = sp.array(counts_new)
                else:
                    counts_anno = sp.array([sp.sum((idx_conf_tcga & is_anno_tcga)[idx_psi_tcga[(mr, ct, dp)]]) if (mr, ct, dp) in idx_psi_tcga else 0 for ct in ctypes])
                    counts_new = sp.array([sp.sum(idx_conf_tcga[idx_psi_tcga[(mr, ct, dp)]]) - sp.sum((idx_conf_tcga & is_anno_tcga)[idx_psi_tcga[(mr, ct, dp)]]) if (mr, ct, dp) in idx_psi_tcga else 0 for ct in ctypes])

                ### plot stats for events by cancer type
                if e == 0:
                    figs['stats_mr%i_dp%.1f' % (mr, dp)] = plt.figure(figsize=(15, 9)) 
                ax = figs['stats_mr%i_dp%.1f' % (mr, dp)].add_subplot(gss['stats'][e / 2, e % 2])
                tickrange1 = sp.arange(ctypes.shape[0])
                ax.bar(tickrange1 + 0.2, counts_anno, 0.6, color=colors, linewidth=0.5)
                ax.bar(tickrange1 + 0.2, counts_new, 0.6, bottom=counts_anno, color=colors, linewidth=0.5, alpha=0.5)
                axs.set_ticks_outer(ax)
                axs.clean_axis(ax)
                ax.set_xticks(tickrange1 + 0.15)
                ax.set_xlim([-1, ctypes.shape[0]])
                ax.set_xticklabels(ctypes, rotation=90, fontsize=10)
                ax.set_title(event_dict[event_type])
                ax.yaxis.grid(True)

                continue ### skip log for now
                ### plot stats for events by histotype (log)
                if e == 0:
                    figs['stats_log_mr%i_dp%.1f' % (mr, dp)] = plt.figure(figsize=(15, 9)) 
                ax = figs['stats_log_mr%i_dp%.1f' % (mr, dp)].add_subplot(gss['stats'][e / 2, e % 2])
                ax.bar(sp.arange(ctypes.shape[0]) + 0.1, sp.log10(counts_anno + 1), 0.8, color=colors, linewidth=0.5)
                ax.bar(sp.arange(ctypes.shape[0]) + 0.1, sp.log10(counts_new + 1), 0.8, bottom=sp.log10(counts_anno + 1), color=colors, linewidth=0.5, alpha=0.3)
                axs.set_ticks_outer(ax)
                axs.clean_axis(ax)
                ax.set_xticks(range(ctypes.shape[0]))
                ax.set_xlim([-1, ctypes.shape[0]])
                ax.set_xticklabels(ctypes, rotation=90, fontsize=10)
                ax.set_title(event_dict[event_type])
                ax.yaxis.grid(True)


    for p in figs:
        print 'plotting %s' % p
        figs[p].tight_layout()
        figs[p].savefig(os.path.join(PLOTDIR, 'event_overview_per_ct_norm_C%i_%s.pdf' % (CONF, p)), format='pdf', bbox_inches='tight')
        figs[p].savefig(os.path.join(PLOTDIR, 'event_overview_per_ct_norm_C%i_%s.png' % (CONF, p)), format='png', bbox_inches='tight')
        plt.close(figs[p])

if __name__ == "__main__":
    main()
