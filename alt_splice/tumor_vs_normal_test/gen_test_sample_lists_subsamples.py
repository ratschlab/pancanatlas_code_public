import sys
import scipy as sp
import os
import numpy.random as npr
npr.seed(23)

t_thresh = 50
n_thresh = 10

sys.path.append('/cluster/home/akahles/git/projects/2013/PanCancerTCGA/rerun2017')
from utils.paths import paths

outdir = os.path.join(paths.basedir, 'rerun2018_alt_splice_difftest_subsampled', 'sample_lists')
if not os.path.exists(outdir):
    os.makedirs(outdir)

whitelist = sp.loadtxt(paths.whitelist, delimiter='\t', dtype='str')
metadata = []
for line in open(paths.metadata, 'r'):
    sl = line.strip().split('\t')
    metadata.append(sl)
metadata = sp.array(metadata, dtype='str')
metaheader = metadata[0, :]
metadata = metadata[1:, :]

### subset to whitelisted samples
aidx = sp.where(metaheader == 'aliquot_id')[0]
tidx = sp.where(metaheader == 'tcga_id')[0]
k_idx = sp.where(sp.in1d(sp.array(['.'.join(x) for x in sp.c_[metadata[:, tidx], metadata[:, aidx]]]), whitelist))[0]
metadata = metadata[k_idx, :]

### get all alignment files
samplelist = sp.loadtxt(os.path.join(paths.basedir_as, 'sample_list.txt'), dtype='str', delimiter='\t')
samples = sp.array([x.split('/')[-1].split('.')[1] for x in samplelist], dtype='str') 

### subset to whitelist
whitelist = sp.loadtxt(paths.whitelist, delimiter='\t', dtype='str')
k_idx = sp.in1d(sp.array([x.split('/')[-2] for x in samplelist]), whitelist)
samples = samples[k_idx]
samplelist = samplelist[k_idx]

for run in range(10):
    outdir_ = os.path.join(outdir, 'run_%i' % (run + 1))
    if not os.path.exists(outdir_):
        os.makedirs(outdir_)

    ### write all tumor sample lists
    cidx = sp.where(metaheader == 'study')[0]
    tidx = sp.where(metaheader == 'is_normal')[0]
    for ctype in sp.unique(metadata[:, cidx]):
        if sp.sum(metadata[:, cidx] == ctype) < 90:
            continue
        ctidx = sp.where((metadata[:, cidx] == ctype) & (metadata[:, tidx] == 'False'))[0]
        if ctidx.shape[0] > t_thresh:
            ctidx = npr.choice(ctidx, t_thresh, replace=False)

        if ctidx.shape[0] == 0:
            continue

        sidx = sp.where(sp.in1d(samples, metadata[ctidx, aidx]))[0]
        sp.savetxt(os.path.join(outdir_, '%s_tumor.txt' % ctype), samplelist[sidx], fmt='%s', delimiter='\t')

    ### write all normal sample lists
    for ctype in sp.unique(metadata[:, cidx]):
        if sp.sum(metadata[:, cidx] == ctype) < 90:
            continue
        ctidx = sp.where((metadata[:, cidx] == ctype) & (metadata[:, tidx] == 'True'))[0]
        if ctidx.shape[0] > n_thresh:
            ctidx = npr.choice(ctidx, n_thresh, replace=False)

        if ctidx.shape[0] == 0:
            continue

        sidx = sp.where(sp.in1d(samples, metadata[ctidx, aidx]))[0]
        sp.savetxt(os.path.join(outdir_, '%s_normal.txt' % ctype), samplelist[sidx], fmt='%s', delimiter='\t')
