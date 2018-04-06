import sys
import os
import h5py
import scipy as sp
import re

sys.path.append('/cluster/home/akahles/git/projects/2013/PanCancerTCGA/rerun2017')
from utils.paths import paths

CONF = 2
event_types = ['exon_skip', 'alt_3prime', 'alt_5prime', 'intron_retention']

### load size factors
libsize = sp.loadtxt(paths.libsize, dtype='str', delimiter='\t')
libsize = libsize[1:, :]
sidx = sp.argsort(libsize[:, 0])
libsize = libsize[sidx, :][:, [0, 1]]
ml = sp.median(libsize[:, 1].astype('float'))
sf = ml / libsize[:, 1].astype('float')


### iterate over events
for et in event_types:
    print 'processing %s events' % et
    IN = h5py.File(os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.counts.hdf5' % (et, CONF)), 'r')
    c = IN['iso1'].chunks
    strains = sp.array([re.sub(r'.aligned$', '', x) for x in IN['strains'][:]])
    sidx = sp.argsort(strains)
    assert sp.all(libsize[:, 0] == strains[sidx])
    rev_sidx = sp.argsort(sidx)

    OUT = h5py.File(os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.min_read_filt.hdf5' % (et, CONF)), 'w')
    OUT.create_dataset(name='strains', data=IN['strains'][:], compression='gzip')

    for ii, i in enumerate(range(0, IN['iso1'].shape[1], c[1])):
        if ii > 0:
            sys.stdout.write('.')
            if ii % 50 == 0:
                sys.stdout.write('%i/%i\n' % (ii, IN['iso1'].shape[1] / c[1]))
            sys.stdout.flush()
        idx = sp.arange(i, min(i+c[1], IN['iso1'].shape[1]))
        tmp1 = IN['iso1'][:, idx][sidx, :]
        tmp2 = IN['iso2'][:, idx][sidx, :]

        tmp3 = sp.c_[(tmp1 * sf[:, sp.newaxis]).max(axis=0), (tmp2 * sf[:, sp.newaxis]).max(axis=0)].min(axis=1)
        for k in [5, 10, 20, 30, 50]:
            if ii == 0:
                OUT.create_dataset(name='minreads_%i' % k, data=(tmp3 >= k), compression='gzip', dtype='bool', maxshape=(None,))
            else:
                ts = OUT['minreads_%i' % k].shape
                OUT['minreads_%i' % k].resize((ts[0] + idx.shape[0],))
                OUT['minreads_%i' % k][ts[0]:] = data=(tmp3 >= k)
        del tmp1, tmp2, tmp3
    OUT.close()
