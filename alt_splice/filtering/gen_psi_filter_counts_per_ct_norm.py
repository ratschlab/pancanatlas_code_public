import sys
import h5py
import scipy as sp
import cPickle
import re

import numpy.random as npr
npr.seed(23)

sys.path.append('/cluster/home/akahles/git/projects/2013/PanCancerTCGA/rerun2017')
import utils.utils as utils
from utils.paths import paths

FILT_THRESH = 40

if len(sys.argv) < 2:
    print >> sys.stderr, 'Usage: %s <counts.hdf5>' % sys.argv[0]
    sys.exit(1)
infname = sys.argv[1]

whitelist = sp.loadtxt(paths.whitelist, delimiter='\t', dtype='str')
whitelist = sp.array([x.split('.')[1] for x in whitelist])

d_psi_t = [0.0, 0.1, 0.3, 0.5]
mr_t = [5, 20, 50]
nan_t = 10

IN = h5py.File(infname, 'r')

### get strain / CT information
(ct_dict, tumor_dict) = utils.get_ct_dict_metatable(paths.metadata, style='pancan_rerun18')
strains = sp.array([x.split('.')[1] for x in IN['strains'][:]])
ctypes = sp.array([ct_dict[x] if x in ct_dict else 'NA' for x in strains])
istumor = sp.array([tumor_dict[x] if x in tumor_dict else 'NA' for x in strains])
ctypes_u, ctypes_cnt = sp.unique(ctypes, return_counts=True)

### reduce to whitelist and only keep tumor samples
fidx = sp.where(~istumor | ~sp.in1d(strains, whitelist))[0]
ctypes[fidx] = 'NA'

### generate a random subset of FILT_THRESH(40) for all ctypes (should be the size of smallest set)
for ct in ctypes_u:
    if ct == 'NA':
        continue
    cidx = sp.where(ctypes == ct)[0]
    if cidx.shape[0] > FILT_THRESH:
        nidx = npr.choice(cidx, cidx.shape[0] - FILT_THRESH, replace=False)
        ctypes[nidx] = 'NA'
ctypes_u, ctypes_cnt = sp.unique(ctypes, return_counts=True)

### load libsize
libsize = sp.loadtxt(paths.libsize, dtype='str', delimiter='\t')
libsize = libsize[1:, :]
libsize[:, 0] = sp.array([x.split('.')[1] for x in libsize[:, 0]])
a, b = sp.where(strains[:, sp.newaxis] == libsize[:, 0])
assert sp.all(strains == libsize[b, 0])
libsize = libsize[b, :]
sfm = sp.median(libsize[:, 1].astype('float'))
sf = sfm / libsize[:, 1].astype('float')

### determine which events to keep
print >> sys.stderr, 'Collecting filter information ...'
chunksize = IN['psi'].chunks[1]

k_idx = dict()
for cc, chunk in enumerate(range(0, IN['psi'].shape[1], chunksize)):

    if cc > 0 and cc % 10 == 0:
        sys.stderr.write('.')
        if cc % 100 == 0:
            sys.stderr.write('%i/%i\n' % (cc, IN['psi'].shape[1] / chunksize + 1))
        sys.stderr.flush()

    tmp_idx = sp.arange(chunk, min(chunk + chunksize, IN['psi'].shape[1]))

    tmp_psi_ = IN['psi'][:, tmp_idx]
    tmp_iso1_ = IN['iso1'][:, tmp_idx] * sf[:, sp.newaxis]
    tmp_iso2_ = IN['iso2'][:, tmp_idx] * sf[:, sp.newaxis]

    for ct in ctypes_u:
        if ct == 'NA':
            continue
        hidx = sp.where((ctypes == ct) & istumor)[0]

        ### no confident events in current chunk
        if hidx.shape[0] < 2:
            continue
        tmp_psi = tmp_psi_[hidx, :] 
        tmp_iso1 = tmp_iso1_[hidx, :]
        tmp_iso2 = tmp_iso2_[hidx, :]
        for mr in mr_t:
            n_idx = sp.c_[tmp_iso1.max(axis=0), tmp_iso2.max(axis=0)].min(axis=1) < mr
            tmp_psi[:, n_idx] = sp.nan
            idx_nn = ~sp.isnan(tmp_psi)
            d_psi = sp.nanmax(tmp_psi, axis=0) - sp.nanmin(tmp_psi, axis=0)
            d_psi[sp.isnan(d_psi)] = 0
            for dp in d_psi_t:
                if not (mr, ct, dp) in k_idx:
                    k_idx[(mr, ct, dp)] = tmp_idx[(sp.sum(idx_nn, axis=0) >= nan_t) & (d_psi >= dp)]
                else:
                    k_idx[(mr, ct, dp)] = sp.r_[k_idx[(mr, ct, dp)], tmp_idx[(sp.sum(idx_nn, axis=0) >= nan_t) & (d_psi >= dp)]]

cPickle.dump((d_psi_t, k_idx), open(re.sub(r'.hdf5$', '', infname) + '.psi_filt_per_ct_normalized.pickle', 'w'), -1)

IN.close()
