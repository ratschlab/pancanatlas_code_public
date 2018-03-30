import sys
import scipy as sp
import scipy.stats as spst
import h5py
import os

sys.path.append('/cluster/home/akahles/git/projects/2013/PanCancerTCGA/rerun2017')
from utils.paths import paths
import utils.samples as samples
import utils.utils as utils

sys.path.append('/cluster/home/akahles/git/projects/2013/PanCancerTCGA/pysrc/utilities')
import names as nm

CONF = 2

use_wl = True
wl_tag = ''
if use_wl:
    wl_tag = '.whitelisted'

outdir = os.path.join(paths.basedir_as, 'outliers')
if not os.path.exists(outdir):
    os.makedirs(outdir)

### get TCGA type dictionary
(ct_dict, is_tumor_dict) = utils.get_ct_dict_metatable(paths.metadata, style='pancan_rerun18')

### walk through all events and check whether a subset of samples consists of extreme outliers
event_types = ['alt_3prime', 'alt_5prime', 'intron_retention', 'exon_skip', 'mutex_exons']

outliers = dict()
out_thresh = 3
out_fact = 10
dpsi = 0.4

lookup = nm.get_lookup_complete()

for et in event_types:
    outlier_cnt = 0
    IN = h5py.File(os.path.join(paths.basedir_as, 'merge_graphs_%s_C%i.counts.hdf5' % (et, CONF)), 'r')
    IN_GT = h5py.File(os.path.join(paths.basedir_as_gtex, 'merge_graphs_%s_C%i.counts.hdf5' % (et, CONF)), 'r')
    strains = sp.array([x.split('.')[1] for x in IN['strains'][:]])

    ### apply whitelist
    if use_wl:
        whitelist = sp.loadtxt(paths.whitelist, delimiter='\t', dtype='str')
        whitelist = sp.array([x.split('.')[1] for x in whitelist])
        widx = sp.in1d(strains, whitelist)
        strains = strains[widx]
    else:
        widx = sp.arange(strains)
        
    ctypes = sp.array([ct_dict[x] for x in strains])
    is_tumor = sp.array([is_tumor_dict[x] for x in strains])
    ctypes_u, ctypes_cnt = sp.unique(ctypes, return_counts=True)
    gene_ids = IN['gene_names'][:]
    gene_names = sp.array([nm.get_ID(x, lookup=lookup) for x in gene_ids])
    gene_idx = IN['gene_idx'][:].astype('int')

    ### only keep cancer types where we have at least 100 samples
    k_idx = sp.where(ctypes_cnt >= 100)[0]
    ctypes_u = ctypes_u[k_idx]

    outliers[et] = dict([(x, []) for x in ctypes_u])
    chunks = IN['psi'].chunks

    for cc, c in enumerate(range(0, IN['psi'].shape[1], chunks[1])): 
        sys.stdout.write('%i/%i - outliers: %i\n' % (cc, IN['psi'].shape[1] / chunks[1], outlier_cnt))
        sys.stdout.flush()
        cidx = sp.arange(c, min(c + chunks[1], IN['psi'].shape[1]))
        tmp = IN['psi'][:, cidx][widx, :]
        tmp_gt = IN_GT['psi'][:, cidx]
        gt_dpsi = sp.nanmax(tmp_gt, axis=0) - sp.nanmin(tmp_gt, axis=0)
        tmp_tn = tmp[~is_tumor, :]
        tn_dpsi = []
        for i in range(tmp_tn.shape[1]):
            ns_idx = sp.argsort(tmp_tn[:, i])
            nnn_idx = ~sp.isnan(tmp_tn[:, i][ns_idx])
            if sp.sum(nnn_idx) > 3:
                tn_dpsi.append(tmp_tn[:, i][ns_idx][nnn_idx][-2] - tmp_tn[:, i][ns_idx][nnn_idx][1])
            else:
                tn_dpsi.append(sp.nanmax(tmp_tn[:, i]) - sp.nanmin(tmp_tn[:, i])) 
        for ct in ctypes_u:
            ctidx = sp.where((ctypes == ct) & is_tumor)[0]
            cnidx = sp.where((ctypes == ct) & ~is_tumor)[0]
            for i in range(tmp.shape[1]):
                nnidx = ~sp.isnan(tmp[ctidx, i])
                if sp.sum(nnidx) < 80:
                    continue

                if not sp.isnan(gt_dpsi[i]) and gt_dpsi[i] > 0.3:
                    continue

                if not sp.isnan(tn_dpsi[i]) and tn_dpsi[i] > 0.3:
                    continue

                ctmp = tmp[ctidx, i][nnidx]
                cdpsi = ctmp.max() - ctmp.min()

                if cdpsi < dpsi:
                    continue

                p25 = spst.scoreatpercentile(ctmp, 25)
                p75 = spst.scoreatpercentile(ctmp, 75)
                iqr = p75 - p25
                cnt = max(sp.sum(ctmp > p75+(out_fact*iqr)), sp.sum(ctmp < p25-(out_fact*iqr)))
                ### outlier too abundant
                if cnt > 100:
                    continue
                ### check relation to normal samples
                if cnidx.shape[0] > 0 and not sp.all(sp.isnan(tmp[cnidx, i])) and sp.nanmin(tmp[cnidx, i]) < ctmp.min():
                    continue

                ### check IQR
                if (iqr > 0) and (cnt >= out_thresh):
                    outliers[et][ct].append([c+i, sp.sum(nnidx), cdpsi, iqr, cnt])
                    outlier_cnt += 1
    IN.close()

    ### assemble output
    header = sp.array(['cancer_type', 'event_type', 'gene_name', 'gene_ids', 'event_id', 'support', 'dpsi', 'iqr', 'outliers'])
    data = []
    for ct in outliers[et]:
        for e in outliers[et][ct]:
            tmp = [ct, et, gene_names[gene_idx[e[0]]], gene_ids[gene_idx[e[0]]], e[0]]
            tmp.extend(e[1:])
            data.append(tmp)
    data = sp.array(data)
    sp.savetxt(os.path.join(outdir, 'outliers_%s%s.tsv' % (et, wl_tag)), sp.r_[header[sp.newaxis, :], data], fmt='%s', delimiter='\t')
