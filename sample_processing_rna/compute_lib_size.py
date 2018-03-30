import sys
import scipy as sp
import h5py
import os
import re

if len(sys.argv) < 2:
    print >> sys.stderr, 'Usage: %s <expression.hdf5>' % sys.argv[0]
    sys.exit(1)
fname = sys.argv[1]

### get list of protein coding genes
coding = sp.loadtxt('/cluster/work/grlab/projects/TCGA/PanCancer/annotation/gencode.v19.annotation.coding_gene_list.tab', dtype='str', delimiter='\t')
coding = coding[:, [0, 4]]
### filter for autosomes
k_idx = sp.where(~sp.in1d(coding[:, 0], sp.array(['chrMT', 'chrX', 'chrY'])))[0]
coding = coding[k_idx, :]
coding = coding[:, 1]

print 'loading expression data from ' + fname
IN = h5py.File(fname, 'r')
expression = IN['counts'][:]
genes = IN['gids'][:]
strains = IN['sids'][:]
IN.close()

k_idx = sp.where(sp.in1d(genes, coding))[0]
genes = genes[k_idx]
expression = expression[k_idx, :]
libsize_uq = sp.array([sp.percentile(x, 75) for x in expression.T])
libsize_tc = expression.sum(axis=0)

s_idx = sp.argsort(libsize_uq)[::-1]

out = open(re.sub('.hdf5$', '', fname) + '.libsize.tsv', 'w')
print >> out, 'sample\tlibsize_75percent\tlibsize_total_count'
for i in s_idx:
    print >> out, '\t'.join([strains[i], str(libsize_uq[i]), str(libsize_tc[i])])
out.close()




