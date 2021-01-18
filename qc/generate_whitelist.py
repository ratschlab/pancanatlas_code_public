import sys
import os
import scipy as sp

FASTQC = # file fastqc_overview_table.pub.tsv.gz available here: https://gdc.cancer.gov/about-data/publications/PanCanAtlas-Splicing-2018
DEGRADATION =  'degradation_scores.tsv.gz'
SEQUENCER = 'Sequencer.tsv.gz'
OUTPUTDIR = # directory for the output to be written to

### we flag a library as low quality, if at least 3 of the following criteria get a 'fail':
###     - basewise_qual_flag
###     - overrep_seq_flag
###     - seqwise_gc_cont_flag
###     - basewise_n_cont_flag
###     - seqwise_qual_flag
### we flag a sample as low quality if the degradation score is larger than Q3 + 1.5xIQR
### we flag a sample as low quality if the GC content is more then 1.5xIQR below Q1 or above Q3
### we flag a sample as low quality if the number of reads is more then 1.5xIQR below Q1 or above Q3
### we exclude a sample if it is flagged for at least three low quality criteria
### we exclude a sample if the degradation score is larger than Q3 + 3xIQR
### we exclude a sample if the GC content is more then 3xIQR below Q1 or above Q3
### we exclude a sample if the number of reads is more then 3xIQR below Q1 or above Q3
### we exclude a sample if it has been sequenced on a different machine than Hiseq 2000


### load and process FASTQC data
data_fq = sp.loadtxt(FASTQC, dtype='str', delimiter='\t')
idx_aid = sp.where(data_fq[0, :] == 'analysis_id')[0][0] 
idx_bqf = sp.where(data_fq[0, :] == 'basewise_qual_flag')[0][0] 
idx_osf = sp.where(data_fq[0, :] == 'overrep_seq_flag')[0][0] 
idx_sgf = sp.where(data_fq[0, :] == 'seqwise_gc_cont_flag')[0][0] 
idx_bnf = sp.where(data_fq[0, :] == 'basewise_n_cont_flag')[0][0]
idx_sqf = sp.where(data_fq[0, :] == 'seqwise_qual_flag')[0][0]
idx_rc = sp.where(data_fq[0, :] == 'readcount')[0][0]
idx_gc = sp.where(data_fq[0, :] == 'seqwise_gc_mean')[0][0]

data_fq = data_fq[1:, :]
data_fq_ids = sp.array([x.split('/')[-1] for x in data_fq[:, idx_aid]])
s_idx = sp.argsort(data_fq_ids)
data_fq_ids = data_fq_ids[s_idx]
data_fq = data_fq[s_idx, :]

data_fq_rc = data_fq[:, idx_rc].astype('int')
data_fq_gc = data_fq[:, idx_gc].astype('float')
data_fq = data_fq[:, [idx_bqf, idx_osf, idx_sgf, idx_bnf, idx_sqf]]
data_fq[data_fq == 'pass'] = '0'
data_fq[data_fq == 'warn'] = '1'
data_fq[data_fq == 'fail'] = '2'
data_fq = data_fq.astype('int')

data_fq_ids, u_cnt = sp.unique(data_fq_ids, return_counts=True)
data_fq_ = sp.zeros((u_cnt.shape[0], data_fq.shape[1]), dtype='int')
data_fq_rc_ = sp.zeros((u_cnt.shape[0], ), dtype='int')
data_fq_gc_ = sp.zeros((u_cnt.shape[0], ), dtype='int')
cum = 0
for i in range(u_cnt.shape[0]):
    if u_cnt[i] == 1:
        data_fq_[i, :] = data_fq[cum, :]
        data_fq_rc_[i] = data_fq_rc[cum]
        data_fq_gc_[i] = data_fq_gc[cum]
    else:
        data_fq_[i, :] = data_fq[cum:cum+u_cnt[i], :].max(axis=0)
        data_fq_rc_[i] = sp.mean(data_fq_rc[cum:cum+u_cnt[i]])
        data_fq_gc_[i] = sp.mean(data_fq_gc[cum:cum+u_cnt[i]])
    cum += u_cnt[i]
data_fq = data_fq_
data_fq_rc = data_fq_rc_
data_fq_gc = data_fq_gc_
f_idx_fq = (sp.sum(data_fq == 2, axis=1) > 2)

q25_rc = sp.percentile(data_fq_rc, 25)
q75_rc = sp.percentile(data_fq_rc, 75)
iqr_rc = q75_rc - q25_rc
x_idx_rc = ((data_fq_rc < q25_rc - 3*iqr_rc) | (data_fq_rc > q75_rc + 3*iqr_rc))
f_idx_rc = ((data_fq_rc < q25_rc - 1.5*iqr_rc) | (data_fq_rc > q75_rc + 1.5*iqr_rc))

q25_gc = sp.percentile(data_fq_gc, 25)
q75_gc = sp.percentile(data_fq_gc, 75)
iqr_gc = q75_gc - q25_gc
x_idx_gc = ((data_fq_gc < q25_gc - 3*iqr_gc) | (data_fq_gc > q75_gc + 3*iqr_gc))
f_idx_gc = ((data_fq_gc < q25_gc - 1.5*iqr_gc) | (data_fq_gc > q75_gc + 1.5*iqr_gc))

### load sequencer data
data_sq = sp.loadtxt(SEQUENCER, dtype='str', delimiter='\t')
data_sq = data_sq[1:, :]
k_idx = (data_sq[:, 3] == 'PAIRED')
data_sq = data_sq[k_idx, :]
_, u_idx = sp.unique(data_sq[:, 2], return_index=True) ### make unique on aliquot 
data_sq = data_sq[u_idx, :]
data_sq = sp.c_[sp.array(['.'.join(x[[0, 2]]) for x in data_sq]), data_sq[:, [1, 3, 4]]]
k_idx = sp.in1d(data_sq[:, 0], data_fq_ids)
data_sq = data_sq[k_idx, :]
s_idx = sp.argsort(data_sq[:, 0])
data_sq = data_sq[s_idx, :]
assert sp.all(data_sq[:, 0] == data_fq_ids)
x_idx_sq = (data_sq[:, 3] != 'Illumina HiSeq 2000')

### load degredation data for subset of the samples
data_dg = sp.loadtxt(DEGRADATION, dtype='str', delimiter='\t')
data_dg[:, 0] = sp.array(['.'.join(x.split('/')[-1].split('.')[:2]) for x in data_dg[:, 0]])
q25_dg = sp.percentile(data_dg[:, 1].astype('float'), 25)
q75_dg = sp.percentile(data_dg[:, 1].astype('float'), 75)
iqr_dg = q75_dg - q25_dg
idx = sp.where(sp.in1d(data_dg[:, 0], data_fq_ids))[0]
assert sp.all(data_dg[idx, 0] == data_fq_ids)
data_dg = data_dg[idx, :]

f_idx_dg = (data_dg[:, 1].astype('float') > q75_dg + 1.5*iqr_dg)
x_idx_dg = (data_dg[:, 1].astype('float') > q75_dg + 3*iqr_dg)

### form final whitelist
k_idx = ((sp.c_[f_idx_fq, f_idx_dg, f_idx_rc, f_idx_gc ].astype('int').sum(axis=1) < 2) & ~x_idx_sq & ~x_idx_dg & ~x_idx_rc & ~x_idx_gc)
sp.savetxt(os.path.join(OUTPUTDIR, 'sample_whitelist.txt'), data_fq_ids[k_idx], fmt='%s')
sp.savetxt(os.path.join(OUTPUTDIR, 'sample_excludelist.txt'), data_fq_ids[~k_idx], fmt='%s')

