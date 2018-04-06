import sys
import scipy as sp
import os
import cPickle

sys.path.append('/cluster/home/akahles/git/projects/2013/PanCancerTCGA/rerun2017')
from utils.paths import *

CONF = 2

events = cPickle.load(open(os.path.join(paths.basedir_as, 'merge_graphs_mutex_exons_C%i.pickle' % CONF), 'r'))

eids = sp.array([':'.join([x.chr, str(x.exons1[0, 0]), str(x.exons1[0, 1]), str(x.exons1[2, 0]), str(x.exons1[2, 1])]) for x in events])
eids_u, eids_i = sp.unique(eids, return_index=True)
_, eids_inv = sp.unique(eids, return_inverse=True)

cPickle.dump(eids_i[eids_inv], open(os.path.join(paths.basedir_as, 'merge_graphs_mutex_exons_C%i.cluster_idx.pickle' % CONF), 'w'), -1)

