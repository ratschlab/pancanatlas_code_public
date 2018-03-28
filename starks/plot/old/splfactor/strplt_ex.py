import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt



import scipy as sp
import h5py

import seaborn as sns

import sys
import pdb

fn_gtex_expression = "/cluster/project/raetsch/lab/03/ICGC/qtl_analysis/dge/gtex_fpkm.tsv.gz"
fn_icgc_expression = "/cluster/project/raetsch/lab/03/ICGC/librarySizeNorm/gene_counts_UCSC_freeze_6_aliquotlevel/whitelist/star_fpkm.tsv.gz"
fn_vars = "/cluster/project/raetsch/lab/03/ICGC/qtl_analysis/hdf5/mergedFiles_ma.hdf5"
fn_indx = "/cluster/project/raetsch/lab/03/ICGC/qtl_analysis/hdf5/midxes.0.1.hdf5"

sns.set_style('white')
def geticgcpeerpheno(fn, gid):
    IN = h5py.File(fn, 'r')
    ix = sp.where(IN['peer_col_header'][:] == gid)[0]
    if 'peer35' in IN.keys():
        peer = IN['peer35'][:,ix]
    else:
        peer = IN['peer15'][:,ix]

    ix = sp.where(IN['col_header'][:] == gid)[0]
    pheno = IN['expression'][:,ix]
    gtids = IN['row_header'][:]
    IN.close()
    return gtids, peer, pheno


def getgtexpheno(gid):
    import gzip
    IN = gzip.open(fn_gtex_expression, 'r')
    header = IN.readline().strip('\n').split('\t')
    for l in IN:
        lSpl = l.strip('\n').split('\t')

        if lSpl[0] == gid:
            IN.close()
            return header,lSpl[1:]
            

def getgenoeinfo(pos):
    chrm = pos.split('_')[0]
    chrm = chrm.strip('chr')
    chrm = int(chrm)
    pos = int(pos.split('_')[1])
    IN = h5py.File(fn_vars,'r')
    allpos = IN['pos'][:]
    ix = sp.where((allpos[:,0] == chrm) & (allpos[:,1] == pos))[0]
    if ix.shape[0] == 0:
        print 'Variant not known'
        sys.exit()
    gt = IN['gt'][ix,:]
    gids = IN['gtid'][:]

    IN.close()
    return gids, gt
    
if __name__ == "__main__":
    sns.set_palette(sns.color_palette("Set2", 10))
    fn_icgc_peer = sys.argv[1]
    gid = sys.argv[2]
    pos = sys.argv[3] ## format ist chr3_49245
    ctype = sys.argv[4]
    import names

    lookup = names.get_lookup_complete()

    gid_trivial = names.get_ID(gid, lookup = lookup)
    peerpheno = geticgcpeerpheno(fn_icgc_peer, gid)

    gtexpheno = getgtexpheno(gid)    
    gtexmap = sp.loadtxt('/cluster/project/raetsch/lab/03/ICGC/gtex_tables/SraRunTable.txt', delimiter = '\t', dtype = 'string', usecols = [6,9])
    gtexmap = dict(gtexmap)

    geno = getgenoeinfo(pos)
    IN = h5py.File(fn_indx, 'r')
    
    rna_gtid = IN['rna_gtid'][:]
    gtidpeer = peerpheno[0][:,1]
    midx_rna = sp.array([sp.where(x == gtidpeer)[0] for x in rna_gtid if x in gtidpeer]).ravel()

    peery = peerpheno[1][midx_rna].ravel()
    fpkmy = peerpheno[2][midx_rna].ravel()
    gt = geno[1].ravel()[IN['idx_dna'][:]]
    tissue = IN['tissue'][:]

    if ctype != 'Pan':
        gt = gt[tissue == ctype]
        tissue = tissue[tissue == ctype]

    iOK = ~sp.isnan(gt)
    peery = peery[iOK]
    fpkmy = fpkmy[iOK]
    gt = gt[iOK]
    tissue = tissue[iOK]

    

    if sp.sum(gt == 2) > sp.sum(gt == 0):
        effsizefpkm = sp.log2(sp.median(fpkmy[gt == gt.max()]) / sp.median(fpkmy[gt == gt.min()]))
    else:
        effsizefpkm = sp.log2(sp.median(fpkmy[gt == gt.min()]) / sp.median(fpkmy[gt == gt.max()]))
    refal = pos.split('_')[2]
    altal = pos.split('_')[3]

    fig = plt.figure(figsize = (20,25)) ### really should be using gridspecs
    gt = gt.astype('int').astype("|S3")
#    gt[gt == '0'] = '%s/%s' % (refal,refal)
#    gt[gt == '1'] = '%s/%s' % (refal,altal)
#    gt[gt == '2'] = '%s/%s' % (altal,altal)

    ax = fig.add_subplot(321)
    ax = sns.stripplot(x=gt,y=peery,hue=tissue, ax = ax,jitter=True)
    ax.set_xticklabels(['%s/%s' % (refal,refal),'%s/%s' % (refal,altal),'%s/%s' % (altal,altal)])
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=4, mode="expand", borderaxespad=0.)
#    ax.set_title('Peer residuals')
#    plt.tight_layout()
    ax = fig.add_subplot(322)
    ax = sns.stripplot(x=gt,y=sp.log10(fpkmy+1),hue=tissue, ax = ax,jitter=True,palette=sns.color_palette("Set2", 20))
    ax = sns.boxplot(x=gt,y=sp.log10(fpkmy+1),ax = ax,fliersize=0)
    ax.set_xticklabels(['%s/%s' % (refal,refal),'%s/%s' % (refal,altal),'%s/%s' % (altal,altal)])
    ax.set_ylabel('log10(FPKM+1)')
    ax.set_ylim(bottom=0.0)
    ax.set_title('FPKM')
    ax.legend_.remove()

    ax = fig.add_subplot(323)
    ax = sns.violinplot(x=gt,y=peery,ax=ax)
    ax.set_xticklabels(['%s/%s' % (refal,refal),'%s/%s' % (refal,altal),'%s/%s' % (altal,altal)])
    ax.set_title('Peer residuals')

    ax = fig.add_subplot(324)
    ax = sns.violinplot(x=gt,y=fpkmy, ax = ax)
    ax.set_xticklabels(['%s/%s' % (refal,refal),'%s/%s' % (refal,altal),'%s/%s' % (altal,altal)])
    ax.set_title('FPKM')


    gtexy = sp.array(gtexpheno[1]).astype('float')    

    gtextissue = sp.array([gtexmap[x] for x in gtexpheno[0][1:]])


    ax = fig.add_subplot(313)

#    ax = sns.stripplot(x=gtextissue, y=gtexy, ax = ax, jitter = True)
    iBrain = sp.array([x.startswith('Brain') for x in gtextissue])
    gtextissue[iBrain] = 'Brain'
    iEso = sp.array([x.startswith('Esophagus') for x in gtextissue])
    gtextissue[iEso] = 'Esophagus'


    ax = sns.stripplot(x= sp.hstack((gtextissue,tissue)), y=sp.hstack((gtexy,fpkmy)), ax = ax, jitter = True,palette=sns.color_palette("Set2", 20))
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
#    plt.tight_layout()
    plt.savefig('%s_%s_%s_%f.pdf' % (gid, gid_trivial, pos, effsizefpkm),bbox_inches='tight')
    IN.close()

