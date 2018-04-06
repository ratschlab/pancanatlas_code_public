import sys
import os
import scipy as sp
import glob

studies = ['BLCA', 'BRCA', 'COAD', 'HNSC', 'KIRC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'PRAD', 'READ', 'STAD', 'THCA', 'UCEC'] 
### directory containing the results of the individual differential tests
BASEDIR = 
event_types = ['alt_3prime', 'alt_5prime', 'exon_skip', 'intron_retention']
event_types = ['mutex_exons']

for event_type in event_types:
    for study in studies:
        print 'processing ' + event_type + ':' + study
        data = []
        for run in range(1, 10):
            data.append(sp.loadtxt(os.path.join(BASEDIR, 'testing_%s-T_vs_%s-N_R%i' % (study, study, run), 'test_results_C2_%s.gene_unique.tsv' % event_type), dtype='str', delimiter='\t'))
        OUTDIR = os.path.join(BASEDIR, 'testing_%s-T_vs_%s-N_MERG' % (study, study))
        if not os.path.exists(OUTDIR):
            os.makedirs(OUTDIR)
        outfname = os.path.join(OUTDIR, 'test_results_C2_%s.gene_unique.tsv' % event_type)

        all_genes = sp.unique([x for d in data for x in d[1:, 1]])
        all_genes_dict = dict([(x, i) for i, x in enumerate(all_genes)])
        
        all_data = sp.ones((all_genes.shape[0], len(data)), dtype='float')

        for r, d in enumerate(data):
            for i in range(1, d.shape[0]):
                all_data[all_genes_dict[d[i, 1]], r] = float(d[i, 3])
        all_data[sp.isnan(all_data)] = 1
        s_idx1 = sp.argsort(all_data, axis=1)
        medians = s_idx1[:, 4]

        data_collect = []
        ### assemble final dataset
        for r, d in enumerate(data):
            for i in range(1, d.shape[0]):
                if medians[all_genes_dict[d[i, 1]]] == r:
                    data_collect.append(d[i, :])

        data_collect = sp.array(data_collect)
        k_idx = sp.where(data_collect[:, 3] != 'nan')[0]
        data_collect = data_collect[k_idx, :]
        s_idx = sp.argsort(data_collect[:, 3].astype('float'))
        data_collect = data_collect[s_idx, :]

        sp.savetxt(outfname, sp.r_[data[0][0, :][sp.newaxis, :], data_collect], fmt='%s', delimiter='\t')
