import sys
import os
import scipy as sp
import cPickle
import h5py

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['font.size'] = 18
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.append('/cluster/home/akahles/git/projects/2013/PanCancerTCGA/colors')
import TCGA_colors as tc 

sys.path.append('/cluster/home/akahles/git/tools/python/viz/')
import axes as axs
sys.path.append('/cluster/home/akahles/git/projects/2013/PanCancerTCGA/rerun2017')
from utils.paths import *

BASEDIR_ANNO = '/cluster/work/grlab/projects/TCGA/PanCancer/rerun_alt_splice_anno'
PLOTDIR = os.path.join(paths.plotdir, 'stats')
if not os.path.exists(PLOTDIR):
    os.makedirs(PLOTDIR)
CONF = 2

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
event_dict = {'exon_skip':'ES',
              'intron_retention':'IR',
              'alt_3prime':'A3',
              'alt_5prime':'A5',
              'mutex_exons':'MX'}


def main():
    figs = dict()
    figs['ctypes'] = plt.figure(figsize=(16, 4))
    figs[''] = plt.figure(figsize=(5, 5))
    gss = dict()
    gss['ctypes'] = gridspec.GridSpec(4, 8) #, wspace=0.0, hspace=0.0)
    gss['stats_log'] = gridspec.GridSpec(2, 2) #, wspace=0.0, hspace=0.0)

    anno = dict()
    is_anno_gtex = dict()
    is_anno_tcga = dict()
    idx_conf_gtex = dict()
    idx_conf_tcga = dict()
    idx_psi_tcga = dict()
    counts_anno = dict()
    counts_new = dict()

    mutex_cluster = None

    for e, event_type in enumerate(event_types):

        print >> sys.stderr, 'Handling %s' % event_type
        
        ### load events detected in annotation only
        anno[event_type] = cPickle.load(open(os.path.join(BASEDIR_ANNO, 'merge_graphs_%s_C3.pickle' % (event_type)), 'r'))
        if isinstance(anno[event_type], tuple):
            anno[event_type] = anno[event_type][0]

        ### load annotation index
        is_anno_gtex[event_type] = cPickle.load(open(os.path.join(paths.basedir_as_gtex, 'merge_graphs_%s_C%i.anno_only.pickle' % (event_type, CONF)), 'r'))
        is_anno_tcga[event_type] = cPickle.load(open(os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.anno_only.pickle' % (event_type, CONF)), 'r'))

        ### load mutex cluster idx
        if event_type == 'mutex_exons':
            mutex_clusters = cPickle.load(open(os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.cluster_idx.pickle' % (event_type, CONF)), 'r'))

        ### load confident events
        IN = h5py.File(os.path.join(paths.basedir_as_gtex, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, CONF)), 'r')
        idx_conf_gtex[event_type] = IN['conf_idx'][:]
        IN.close()
        IN = h5py.File(os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.counts.hdf5' % (event_type, CONF)), 'r')
        idx_conf_tcga[event_type] = sp.zeros((is_anno_tcga[event_type].shape[0],), dtype='bool')
        idx_conf_tcga[event_type][IN['conf_idx'][:]] = 1
        IN.close()

        ### load psi filtered events
        idx_psi_tcga = cPickle.load(open(os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.counts.psi_filt_per_ct_normalized.pickle' % (event_type, CONF)), 'r'))[1]

        ### get all histotypes
        ctypes = sp.unique([x[1] for x in idx_psi_tcga.keys()])
        colors = tc.get_color_scheme('study', labels=ctypes)

        ### get counts
        dp = 0.3
        mr = 20
        if event_type == 'mutex_exons':
            counts_anno_ = []
            counts_new_ = []
            for ct in ctypes:
                if not (mr, ct, dp) in idx_psi_tcga:
                    counts_anno_.append(0)
                    counts_new_.append(0)
                    continue
                tmp_idx = sp.zeros((is_anno_tcga[event_type].shape[0],), dtype='bool')
                tmp_idx[idx_psi_tcga[(mr, ct, dp)]] = 1
                counts_anno_.append(sp.unique(mutex_clusters[idx_conf_tcga[event_type] & is_anno_tcga[event_type] & tmp_idx]).shape[0])
                counts_new_.append(sp.unique(mutex_clusters[idx_conf_tcga[event_type] & ~is_anno_tcga[event_type] & tmp_idx]).shape[0])
            counts_anno[event_type] = sp.array(counts_anno_)
            counts_new[event_type] = sp.array(counts_new_)
        else:
            counts_anno[event_type] = sp.array([sp.sum((idx_conf_tcga[event_type] & is_anno_tcga[event_type])[idx_psi_tcga[(mr, ct, dp)]]) if (mr, ct, dp) in idx_psi_tcga else 0 for ct in ctypes])
            counts_new[event_type] = sp.array([sp.sum(idx_conf_tcga[event_type][idx_psi_tcga[(mr, ct, dp)]]) - sp.sum((idx_conf_tcga[event_type] & is_anno_tcga[event_type])[idx_psi_tcga[(mr, ct, dp)]]) if (mr, ct, dp) in idx_psi_tcga else 0 for ct in ctypes])

    colors = ['#d7191c', '#fdae61', '#abd9e9', '#2c7bb6', '#999999']

    ### plot stats for event types by cancer type
    fig = plt.figure(figsize=(16, 8))
    gs = gridspec.GridSpec(4, 8)
    for c,ct in enumerate(ctypes):
        data = []
        for et in event_types:
            data.append(counts_anno[et][c])
        ax = fig.add_subplot(gs[c / 8, c % 8])
        ax.pie(data, radius=1, colors=colors, wedgeprops={'alpha':0.7, 'linewidth':0.1}, autopct=lambda(p): '{:.1f}%'.format(p)) 
        ax.set_title(ct)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTDIR, 'event_type_overview_per_ct_C%i.pdf' % (CONF)), format='pdf', bbox_inches='tight')
    plt.savefig(os.path.join(PLOTDIR, 'event_type_overview_per_ct_C%i.png' % (CONF)), format='png', bbox_inches='tight')
    plt.close(fig)

    ### plot stats for annotated events
    fig = plt.figure(figsize=(5, 5))
    data = []
    for et in event_types:
        data.append(anno[et].shape[0])
    ax = fig.add_subplot(111)
    labels = [event_dict[_] for _ in event_types]
    ax.pie(data, labels=labels, radius=1, labeldistance=1.10, colors=colors, wedgeprops={'alpha':0.7, 'linewidth':0.1}, autopct=lambda(p): '{:.1f}%'.format(p)) 
    ax.set_title('Events in Annotation')
    plt.savefig(os.path.join(PLOTDIR, 'event_type_overview_anno_C%i.pdf' % (CONF)), format='pdf', bbox_inches='tight')
    plt.savefig(os.path.join(PLOTDIR, 'event_type_overview_anno_C%i.png' % (CONF)), format='png', bbox_inches='tight')
    plt.close(fig)

    ### plot stats for gtex events
    fig = plt.figure(figsize=(5, 5))
    data = []
    for et in event_types:
        data.append(idx_conf_gtex[et].shape[0])
    ax = fig.add_subplot(111)
    labels = [event_dict[_] for _ in event_types]
    ax.pie(data, labels=labels, radius=1, labeldistance=1.10, colors=colors, wedgeprops={'alpha':0.7, 'linewidth':0.1}, autopct=lambda(p): '{:.1f}%'.format(p)) 
    ax.set_title('Events in GTEx')
    plt.savefig(os.path.join(PLOTDIR, 'event_type_overview_gtex_C%i.pdf' % (CONF)), format='pdf', bbox_inches='tight')
    plt.savefig(os.path.join(PLOTDIR, 'event_type_overview_gtex_C%i.png' % (CONF)), format='png', bbox_inches='tight')
    plt.close(fig)

    ### plot stats for tcga events
    fig = plt.figure(figsize=(5, 5))
    data = []
    for et in event_types:
        data.append(sp.sum(idx_conf_tcga[et]))
    ax = fig.add_subplot(111)
    labels = [event_dict[_] for _ in event_types]
    ax.pie(data, labels=labels, radius=1, labeldistance=1.10, colors=colors, wedgeprops={'alpha':0.7, 'linewidth':0.1}, autopct=lambda(p): '{:.1f}%'.format(p)) 
    ax.set_title('Events in TCGA')
    plt.savefig(os.path.join(PLOTDIR, 'event_type_overview_tcga_C%i.pdf' % (CONF)), format='pdf', bbox_inches='tight')
    plt.savefig(os.path.join(PLOTDIR, 'event_type_overview_tcga_C%i.png' % (CONF)), format='png', bbox_inches='tight')
    plt.close(fig)

    ### plot stats for all filtered tcga events
    fig = plt.figure(figsize=(5, 5))
    data = []
    for et in event_types:
        ### load psi filtered events
        idx_psi_tcga = cPickle.load(open(os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.counts.psi_filt.pickle' % (et, CONF)), 'r'))[1]
        if et == 'mutex_exons':
            data.append(sp.unique(mutex_clusters[idx_psi_tcga[(20, 'all', 0.0)]]).shape[0])
        else:
            data.append(idx_psi_tcga[(20, 'all', 0.0)].shape[0])
    ax = fig.add_subplot(111)
    labels = [event_dict[_] for _ in event_types]
    ax.pie(data, labels=labels, radius=1, labeldistance=1.10, colors=colors, wedgeprops={'alpha':0.7, 'linewidth':0.1}, autopct=lambda(p): '{:.1f}%'.format(p)) 
    ax.set_title('Events in TCGA - Filtered')
    plt.savefig(os.path.join(PLOTDIR, 'event_type_overview_tcga_filt_C%i.pdf' % (CONF)), format='pdf', bbox_inches='tight')
    plt.savefig(os.path.join(PLOTDIR, 'event_type_overview_tcga_filt_C%i.png' % (CONF)), format='png', bbox_inches='tight')
    plt.close(fig)

    ### plot stats for all filtered tcga events that are annotated
    fig = plt.figure(figsize=(5, 5))
    data = []
    for et in event_types:
        ### load psi filtered events
        idx_psi_tcga = cPickle.load(open(os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.counts.psi_filt.pickle' % (et, CONF)), 'r'))[1]
        if et == 'mutex_exons':
            ti = sp.zeros((is_anno_tcga[et].shape[0], ), dtype='bool')
            ti[idx_psi_tcga[(5, 'all', 0.0)]] = 1
            data.append(sp.unique(mutex_clusters[is_anno_tcga[et] & ti]).shape[0])
        else:
            data.append(sp.sum(is_anno_tcga[et][idx_psi_tcga[(5, 'all', 0.0)]]))
    ax = fig.add_subplot(111)
    labels = [event_dict[_] for _ in event_types]
    ax.pie(data, labels=labels, radius=1, labeldistance=1.10, colors=colors, wedgeprops={'alpha':0.7, 'linewidth':0.1}, autopct=lambda(p): '{:.1f}%'.format(p)) 
    ax.set_title('Events in TCGA - Filtered (Anno)')
    plt.savefig(os.path.join(PLOTDIR, 'event_type_overview_tcga_filt_anno_C%i.pdf' % (CONF)), format='pdf', bbox_inches='tight')
    plt.savefig(os.path.join(PLOTDIR, 'event_type_overview_tcga_filt_anno_C%i.png' % (CONF)), format='png', bbox_inches='tight')
    plt.close(fig)


    ### ployt stats for all filtered tcga events that are not annotated
    fig = plt.figure(figsize=(5, 5))
    data = []
    for et in event_types:
        ### load psi filtered events
        idx_psi_tcga = cPickle.load(open(os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.counts.psi_filt.pickle' % (et, CONF)), 'r'))[1]
        if et == 'mutex_exons':
            ti = sp.zeros((is_anno_tcga[et].shape[0], ), dtype='bool')
            ti[idx_psi_tcga[(20, 'all', 0.0)]] = 1
            data.append(sp.unique(mutex_clusters[~is_anno_tcga[et] & ti]).shape[0])
        else:
            data.append(sp.sum(~is_anno_tcga[et][idx_psi_tcga[(20, 'all', 0.0)]]))
    ax = fig.add_subplot(111)
    labels = [event_dict[_] for _ in event_types]
    ax.pie(data, labels=labels, radius=1, labeldistance=1.10, colors=colors, wedgeprops={'alpha':0.7, 'linewidth':0.1}, autopct=lambda(p): '{:.1f}%'.format(p)) 
    ax.set_title('Events in TCGA - Filtered (Novel)')
    plt.savefig(os.path.join(PLOTDIR, 'event_type_overview_tcga_filt_novel_C%i.pdf' % (CONF)), format='pdf', bbox_inches='tight')
    plt.savefig(os.path.join(PLOTDIR, 'event_type_overview_tcga_filt_novel_C%i.png' % (CONF)), format='png', bbox_inches='tight')
    plt.close(fig)

    

if __name__ == "__main__":
    main()
