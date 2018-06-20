import sys
import pdb
import pysam
import time
import re
import scipy as sp
import h5py
import cPickle
import os

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-a', '--annotation', dest='anno', metavar='FILE', help='annotation file in GTF/GFF3 format', default='-')
    required.add_option('-o', '--outfile', dest='outfile', metavar='FILE', help='outfile to store counts in tab delimited format [stdin]', default='-')
    required.add_option('-A', '--alignment', dest='alignment', metavar='FILE', help='alignment in sam or bam format [stdin - sam]', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-F', '--fields', dest='fields', metavar='STRING', help='annotation fields [exon], comma separated', default='exon')
    optional.add_option('-f', '--filters', dest='filters', metavar='STRING', help='file containing filter maps in hdf5 [-]', default='-')
    optional.add_option('-n', '--filternames', dest='filternames', metavar='STRING', help='list of filter names to use, comma separated, names must be present in t filter hdf5 [names in hdf5 in lex order]', default='-')
    optional.add_option('-t', '--filtertypes', dest='filtertypes', metavar='STRING', help='list of filter types to use, comma separated, either one or same number as filternames, possible types: any, start, all [any]', default='-')
    optional.add_option('-c', '--filtercombs', dest='filtercombs', metavar='STRING', help='list of filter-index combinations: 0,2,4:0,1:... (index relative to filter name list) [one filter in hdf5 at a time]', default='-')
    optional.add_option('-m', '--mask_gene_overlap', dest='mask_gene_overlap', action='store_true', help='mask genomic positions that are annotated with different genes [off]', default=False)
    optional.add_option('-M', '--mask_alternative_overlap', dest='mask_alternative_overlap', action='store_true', help='mask genomic positions that are annotated with both intronic and exonic positions [off]', default=False)
    optional.add_option('-b', '--bam_force', dest='bam_force', action='store_true', help='force BAM as input even if file ending is different from .bam - does not work for STDIN', default=False)
    optional.add_option('-B', '--best_only', dest='best_only', action='store_true', help='count only the best alignment per read [off]', default=False)
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    
    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    return options

def compress_g(g, idx2gene):
    """Find reduced g"""

    g = sorted(g, key = lambda x: len(idx2gene[x]))[::-1]
    g_ = [g[0]]
    seen = idx2gene[g[0]]

    for gg in g[1:]:
        if not all([i in seen for i in idx2gene[gg]]):
            g_.append(gg)
            seen += idx2gene[gg]
    
    return sp.array(g_)

def valid_after_filter(filtermap, filtertype, positions):
    """Description"""

    if filtertype == 'all':
        return not sp.all(filtermap[:, positions])
    elif filtertype == 'start':
        return not filtermap[:, positions[0]]
    elif filtertype == 'any':
        return not sp.any(filtermap[:, positions])
    else:
        return False

def get_filter_settings(options):
    """Parse filter settings from command line options.""" 
    
    if options.filternames != '-':
        filter_names = options.filternames.split(',')
    else:
        hdf_in = h5py.File(options.filters, 'r')
        filter_names = sorted(hdf_in.keys())
        hdf_in.close()

    if options.filtercombs != '-':
        filter_combs = []
        for fc in options.filtercombs.split(':'):
            filter_combs.append(fc.split(','))
            filter_combs[-1] = [int(x) for x in filter_combs[-1]]
    else:
        filter_combs = [[x] for x in range(len(filter_names))]

    if options.filtertypes == '-':
        filter_types = ['any'] * len(filter_names)
    else:
        ft = options.filtertypes.split(',')
        if len(ft) == 1:
            filter_types = [ft[0]] * len(filter_names)
        else:
            assert(len(ft) == len(filter_names))
            filter_types = ft
    
    return (filter_names, filter_combs, filter_types)
    

def get_tags_gff(tagline):
    """Extract tags from given tagline"""

    tags = dict()
    for t in tagline.split(';'):
        tt = t.split('=')
        tags[tt[0]] = tt[1]
    return tags


def get_tags_gtf(tagline):
    """Extract tags from given tagline"""

    tags = dict()
    for t in tagline.strip(';').split(';'):
        tt = t.strip(' ').split(' ')
        tags[tt[0]] = tt[1].strip('"')
    return tags


def parse_anno_from_gff3(options, contigs):
    """This function reads the gff3 input file and returns the information in an
       internal data structure"""

    anno = dict()
    idx2gene = dict()
    gene2idx = dict()

    if options.verbose:
        print >> sys.stderr, "Parsing annotation from %s ..." % options.anno
    
    ### initial run to get the transcript to gene mapping
    if options.verbose:
        print >> sys.stderr, "... init structure"

    trans2gene = dict() ### dict with: keys = transcript IDs, values = gene IDs
    for line in open(options.anno, 'r'):
        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        if sl[2] in ['mRNA', 'transcript', 'mrna', 'miRNA', 'tRNA', 'snRNA', 'snoRNA', 'ncRNA', 'mRNA_TE_gene', 'rRNA', 'pseudogenic_transcript', 'transposon_fragment']:
            tags = get_tags_gff(sl[8])
            trans2gene[tags['ID']] = tags['Parent']

    ### init genome structure
    for c in contigs:
        if options.verbose:
            print >> sys.stderr, 'reserving memory for contig %s of len %s' % (c, contigs[c])
        anno[c] = sp.zeros((contigs[c] + 1,), dtype = 'int32')

    ### init list of  considered GFF fields
    fields = options.fields.split(',')

    ### generate a list of exons with attached gene/transcript information
    ### one list per chromsome
    counter = 1
    gene_counter = 2 ### 0 is default for no coverage and 1 is mask for overlap

    exons = dict() # contains the exon list per transcript, only need this for mask_alternative_overlap

    t0 = time.time()
    for line in open(options.anno, 'r'):
        if options.verbose and counter % 10000 == 0:
            print >> sys.stderr, '.',
            if counter % 100000 == 0:
                t1 = time.time() - t0
                print >> sys.stderr, "%i - took %.2f secs" % (counter, t1)
                t0 = time.time()
        counter += 1        

        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        
        if not sl[2] in fields:
            continue

        tags = get_tags_gff(sl[8])
        if sl[2] == 'exon':
            trans_id = tags['Parent']
            gene_id = trans2gene[trans_id]
        else:
            print >> sys.stderr, 'Currently only >exon< is supported'
            sys.exit(1)

        if not gene2idx.has_key(tuple([gene_id])):
            gene2idx[tuple([gene_id])] = gene_counter
            idx2gene[gene_counter] = tuple([gene_id])
            gene_counter += 1

        ### store for each position of the transcriptome a tuple containing all overlapping gene IDs
        ### assume positions are 1 based and in closed intervals
        try:
            start = int(sl[3]) - 1
        except ValueError:
            start = 0
        try:
            stop = int(sl[4])
        except ValueError:
            stop = 1

        if not sl[0] in exons:
            exons[sl[0]] = dict()

        if options.mask_alternative_overlap:
            try:
                exons[sl[0]][trans_id].append([start, stop])
            except KeyError:
                exons[sl[0]][trans_id] = [[start, stop]]

        ### check, if there is already a different gene ID present, form a combination ID
        if sp.any(anno[sl[0]][start:stop] > 0):
            for p in range(start, stop):
                if anno[sl[0]][p] == 0:
                    new_set = tuple([gene_id])
                else:
                    new_set = tuple(set(idx2gene[anno[sl[0]][p]]) | set([gene_id]))
                try:
                    anno[sl[0]][p] = gene2idx[new_set]
                except KeyError:
                    anno[sl[0]][p] = gene_counter
                    gene2idx[new_set] = gene_counter
                    idx2gene[gene_counter] = new_set
                    gene_counter += 1
        else:
            anno[sl[0]][start:stop] = sp.array([gene2idx[tuple([gene_id])]] * (stop - start), dtype = 'int32')
    if options.verbose:
        print >> sys.stderr, "... done"

    ### mask all positions in the genome, where we have more than one annotated gene
    if options.mask_gene_overlap:
        total_pos = 0
        total_masked = 0
        if options.verbose:
            print >> sys.stderr, '\nMasking positions due to gene overlap:'
        for c in anno:
            masked_pos = 0
            p_idx = sp.where(anno[c] > 1)[0]
            pos = p_idx.shape[0]
            for p in p_idx:
                if len(idx2gene[anno[c][p]]) > 1:
                    anno[c][p] = 1
                    masked_pos += 1
            total_pos += pos
            total_masked += masked_pos
            if options.verbose:
                print >> sys.stderr, '\t%s: %i (%i) masked (total) - %.2f %%' % (c, masked_pos, pos, masked_pos / float(max(1, pos)) * 100)
        if options.verbose:
            print >> sys.stderr, "Total positions: %i\nMasked positions: %i (%.2f %%)" % (total_pos, total_masked, total_masked / float(max(1, total_pos)) * 100)
            print >> sys.stderr, "... done"

    ### mask all positions in the genome, where exonic and intronic positions are annotated
    if options.mask_alternative_overlap:
        if options.verbose:
            print >> sys.stderr, '\nMasking positions due to exon/intron overlap:'
        for c in exons:
            masked_pos = 0
            for t in exons[c]:
                if len(exons[c][t]) < 2:
                    continue
                ### pre-process exon
                tmp = sp.array(exons[c][t], dtype='int')
                s_idx = sp.argsort(tmp[:, 0])
                tmp = tmp[s_idx, :]
                ### mask positions that are intronic and exonic
                for e in range(1, tmp.shape[0]):
                    p_idx = sp.where(anno[c][tmp[e - 1, 1] + 1:tmp[e, 0]] > 1)[0]
                    if p_idx.shape[0] > 0:
                        anno[c][p_idx + tmp[e - 1, 1] + 1] = 1
                        masked_pos += p_idx.shape[0]
            total_masked += masked_pos
            if options.verbose:
                print >> sys.stderr, '\t%s: %i pos masked' % (c, masked_pos)
        if options.verbose:
            print >> sys.stderr, 'Masked positions: %i' % total_masked
            print >> sys.stderr, "... done"

    
    if options.verbose:
        print >> sys.stderr, "Storing exon array in HDF5 %s ..." % (options.anno_hdf5 + '.exons.hdf5')

    ### store annotation in hdf5
    hdf_out = h5py.File(options.anno_hdf5 + '.exons.hdf5', 'w')
    for c in anno.keys():
        hdf_out.create_dataset(name = c, data = anno[c])
    hdf_out.close()

    if options.verbose:
        print >> sys.stderr, "... pickling gene ID map"

    cPickle.dump((idx2gene, gene2idx), open(options.anno_hdf5 + '.pickle', 'w'))

    if options.verbose:
        print >> sys.stderr, "... done"

    return (anno, idx2gene, gene2idx)

def parse_anno_from_gtf(options, contigs):
    """This function reads the gtf input file and returns the information in an
       internal data structure"""

    anno = dict()
    idx2gene = dict()
    gene2idx = dict()

    if options.verbose:
        print >> sys.stderr, "Parsing annotation from %s ..." % options.anno
    
    ### init genome structure
    for c in contigs:
        if options.verbose:
            print >> sys.stderr, 'reserving memory for chr %s of len %s' % (c, contigs[c])
        anno[c] = sp.zeros((contigs[c] + 1, ), dtype = 'int32')

    ### init list of  considered GFF fields
    fields = options.fields.split(',')

    ### generate a list of exons with attached gene/transcript information
    ### one list per chromsome
    counter = 1
    gene_counter = 2 ### 0 is default for no coverage and 1 is mask for overlap

    exons = dict()

    t0 = time.time()
    for line in open(options.anno, 'r'):
        if options.verbose and counter % 10000 == 0:
            print >> sys.stderr, '.',
            if counter % 100000 == 0:
                t1 = time.time() - t0
                print >> sys.stderr, "%i - took %.2f secs" % (counter, t1)
                t0 = time.time()
        counter += 1        

        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        
        if not sl[2] in fields:
            continue

        if sl[2] != 'exon':
            print >> sys.stderr, 'Currently only >exon< is supported'
            sys.exit(1)

        tags = get_tags_gtf(sl[8])
        gene_id = tags['gene_id']
        trans_id = tags['transcript_id']

        if not gene2idx.has_key(tuple([gene_id])):
            gene2idx[tuple([gene_id])] = gene_counter
            idx2gene[gene_counter] = tuple([gene_id])
            gene_counter += 1

        try:
            start = int(sl[3]) - 1
        except ValueError:
            start = 0
        try:
            stop = int(sl[4])
        except ValueError:
            stop = 1

        chrm = sl[0]
        if chrm == 'chrM_rCRS':
            chrm = 'chrM'

        if not chrm in exons:
            exons[chrm] = dict()

        if options.mask_alternative_overlap:
            try:
                exons[chrm][trans_id].append([start, stop])
            except KeyError:
                exons[chrm][trans_id] = [[start, stop]]

        ### check, if there is already a different gene ID present, form a combination ID
        if sp.any(anno[chrm][start:stop] > 0):
            for p in range(start, stop):
                if anno[chrm][p] == 0:
                    new_set = tuple([gene_id])
                else:
                    new_set = tuple(set(idx2gene[anno[chrm][p]]) | set([gene_id]))
                try:
                    anno[chrm][p] = gene2idx[new_set]
                except KeyError:
                    anno[chrm][p] = gene_counter
                    gene2idx[new_set] = gene_counter
                    idx2gene[gene_counter] = new_set
                    gene_counter += 1
        else:
            anno[chrm][start:stop] = sp.array([gene2idx[tuple([gene_id])]] * (stop - start), dtype = 'int32')
    if options.verbose:
        print >> sys.stderr, "... done"

    ### mask all positions in the genome, where we have more than one annotated gene
    if options.mask_gene_overlap:
        total_pos = 0
        total_masked = 0
        if options.verbose:
            print >> sys.stderr, '\nMasking positions due to gene overlap:'
        for c in anno:
            masked_pos = 0
            p_idx = sp.where(anno[c] > 1)[0]
            pos = p_idx.shape[0]
            #print >> sys.stderr, 'found %i positions' % p_idx.shape[0]
            for p in p_idx:
                if len(idx2gene[anno[c][p]]) > 1:
                    anno[c][p] = 1
                    masked_pos += 1
            total_pos += pos
            total_masked += masked_pos
            if options.verbose:
                print >> sys.stderr, '\t%s: %i (%i) masked (total) - %.2f %%' % (c, masked_pos, pos, masked_pos / float(max(1, pos)) * 100)
        if options.verbose:
            print >> sys.stderr, "Total positions: %i\nMasked positions: %i (%.2f %%)" % (total_pos, total_masked, total_masked / float(max(1, total_pos)) * 100)
            print >> sys.stderr, "... done"

    ### mask all positions in the genome, where exonic and intronic positions are annotated
    if options.mask_alternative_overlap:
        total_masked = 0
        if options.verbose:
            print >> sys.stderr, '\nMasking positions due to exon/intron overlap:'
        for c in exons:
            masked_pos = 0
            for t in exons[c]:
                if len(exons[c][t]) < 2:
                    continue
                ### pre-process exon
                tmp = sp.array(exons[c][t], dtype='int')
                s_idx = sp.argsort(tmp[:, 0])
                tmp = tmp[s_idx, :]
                ### mask positions that are intronic and exonic
                for e in range(1, tmp.shape[0]):
                    p_idx = sp.where(anno[c][tmp[e - 1, 1] + 1:tmp[e, 0]] > 1)[0]
                    if p_idx.shape[0] > 0:
                        anno[c][p_idx + tmp[e - 1, 1] + 1] = 1
                        masked_pos += p_idx.shape[0]
            total_masked += masked_pos
            if options.verbose:
                print >> sys.stderr, '\t%s: %i pos masked' % (c, masked_pos)
        if options.verbose:
            print >> sys.stderr, 'Masked positions: %i' % total_masked
            print >> sys.stderr, "... done"

    if options.verbose:
        print >> sys.stderr, "Storing exon array in HDF5 %s ..." % (options.anno_hdf5 + '.exons.hdf5')

    ### store annotation in hdf5
    hdf_out = h5py.File(options.anno_hdf5 + '.exons.hdf5', 'w')
    for c in anno.keys():
        hdf_out.create_dataset(name = c, data = anno[c])
    hdf_out.close()

    if options.verbose:
        print >> sys.stderr, "... pickling gene ID map"

    cPickle.dump((idx2gene, gene2idx), open(options.anno_hdf5 + '.pickle', 'w'))

    if options.verbose:
        print >> sys.stderr, "... done"

    return (anno, idx2gene, gene2idx)

def read_header(options, infile):
    """Parses the alignment header and extracts contig information"""

    contigs = dict()
    line = ''
    if options.is_bam:
        #chrm = infile.getrname(line.tid).replace('chr', '')
        for i in range(len(infile.references)):
            if infile.references[i] == 'chrM_rCRS':
                chr_key = 'chrM'
            else:
                chr_key = infile.references[i]

            if contigs.has_key(chr_key):
                if not contigs[chr_key] == infile.lengths[i]:
                    print >> sys.stderr, "Headers in BAM files have inconsistent contig lengths. Stopping ..."
                    sys.exit(1)
            else:
                contigs[chr_key] = infile.lengths[i]
    else:
        for line in infile:
            if not line[0] == '@':
                if len(contigs) == 0:
                    print >> sys.stderr, "No header found in %s. Stopping." % file
                    sys.exit(1)
                else:
                    break

            sl = line.strip().split('\t')

            if not sl[0] == '@SQ':
                continue

            if sl[1][3:] == 'chrM_rCRS':
                chr_key = 'chrM'
            else:
                chr_key = sl[1][3:]
            if contigs.has_key(chr_key):
                if not contigs[chr_key] == int(sl[2][3:]):
                    print >> sys.stderr, "Headers in BAM files have inconsistent contig lengths. Stopping ..."
                    sys.exit(1)
            else:
                contigs[chr_key] = int(sl[2][3:])
                        
    return (contigs, line)

def compress_counts(count_list, genes):
    """Takes a list of gene IDs and compresses them to a list of tuples"""

    a = 0
    g = 0
    compressed_list = []

    print >> sys.stderr, " [compressing gene list] ",

    while g < len(genes):
        while g < len(genes) and (a == len(count_list) or genes[g] < count_list[a]):
            g += 1
        if g < len(genes):
            b = a
            while a < len(count_list) and genes[g] == count_list[a]:
                a += 1
            compressed_list.append([genes[g], a - b])
            g += 1
    return compressed_list
    
def condense_compressed_counts(compressed_counts):

    t0 = time.time()
    for idx in range(len(compressed_counts)):
        compressed_counts[idx] = sorted(compressed_counts[idx], key = lambda x: x[0])
        for i in range(1, len(compressed_counts[idx])):
            if compressed_counts[idx][i-1][0] == compressed_counts[idx][i][0]:
                compressed_counts[idx][i][1] += compressed_counts[idx][i-1][1]
                compressed_counts[idx][i-1][1] = -1
        compressed_counts[idx] = [x for x in compressed_counts[idx] if x[1] >= 0]
        t1 = time.time() - t0
        print >> sys.stderr, "... done. took %.2f secs" % t1

    return compressed_counts


def main():
    """Main Program Procedure"""
    
    options = parse_options(sys.argv)
    contigs = dict()

    options.anno_hdf5 = options.anno
    if options.mask_gene_overlap:
        options.anno_hdf5 += '.mask_go'
    if options.mask_alternative_overlap:
        options.anno_hdf5 += '.mask_ao'

    time_total = time.time()

    ### get filters
    filters = []
    if options.filters != '-':
        ### get filter names
        (filter_names, filter_combs, filter_types) = get_filter_settings(options)

        ### subset to filter names that occur in filter_combs
        filter_combs_flat = list(set([j for sublist in filter_combs for j in sublist]))
        filter_names = [filter_names[i] for i in range(len(filter_names)) if i in filter_combs_flat]
        filter_types = [filter_types[i] for i in range(len(filter_types)) if i in filter_combs_flat]
        filter_combs = [[filter_combs_flat.index(x) for x in j] for j in filter_combs]

        hdf_in = h5py.File(options.filters, 'r')
        for fn in filter_names:
            filters.append(dict())
            for c in hdf_in[fn]:
                filters[-1][c] = hdf_in[fn][c][:]
        hdf_in.close()
    else:
        filter_names = []
        filter_combs = []
        filter_types = []

    ### iterate over alignment file(s)
    for fname in options.alignment.split(','):
        options.is_bam = False
        ### open file stream
        if fname == '-':
            infile = sys.stdin
        elif (len(fname) > 3 and fname[-3:] == 'bam') or options.bam_force:
            infile = pysam.Samfile(fname, 'rb')
            options.is_bam = True
        else:
            infile = open(fname, 'r')

        if options.verbose:
            if options.alignment == '-':
                print >> sys.stderr, "Reading alignment from stdin\n"
            else:
                print >> sys.stderr, "Reading alignment from %s\n" % options.alignment

        ### get contigs from alignment data
        if len(contigs) == 0:
            (contigs, lastline) = read_header(options, infile)
            ### TODO handle lastline (line after header) for SAM input

            ### parse annotation into memory or read from hdf5 is availabl
            if os.path.isfile(options.anno_hdf5 + '.pickle') and os.path.isfile(options.anno_hdf5 + '.exons.hdf5'):
                if options.verbose:
                    t0 = time.time()
                    print >> sys.stderr, 'Loading annotation from %s ...' % (options.anno_hdf5 + '.pickle')
                (idx2gene, gene2idx) = cPickle.load(open(options.anno_hdf5 + '.pickle', 'r'))

                anno = dict()
                hdf_in = h5py.File(options.anno_hdf5 + '.exons.hdf5', 'r')
                for c in hdf_in:
                    anno[c] = hdf_in[c][:]
                    if options.verbose:
                        t1 = time.time() - t0
                        print >> sys.stderr, "... %s took %i secs" % (c, t1)
                        t0 = time.time()
                hdf_in.close()
            else:
                if options.anno[-4:] == 'gff3' or options.anno[-3:] == 'gff':
                    ### read annotation from GFF3
                    (anno, idx2gene, gene2idx) = parse_anno_from_gff3(options, contigs)
                else:
                    ### read annotation from GTF
                    (anno, idx2gene, gene2idx) = parse_anno_from_gtf(options, contigs)

        ### count reads
        counter = 1
        t0 = time.time()
        tmp_count = [[] for i in range(1 + len(filter_combs))]
        compressed_counts = [[] for i in range(1 + len(filter_combs))]
        genes = sorted(idx2gene.keys())
        for line in infile:
            if counter % 10000 == 0:
                print >> sys.stderr, '.',
                if counter % 100000 == 0:
                    if len(tmp_count[0]) > 5000000:
                        for idx in range(len(tmp_count)):
                            compressed_counts[idx].extend(compress_counts(sorted(tmp_count[idx]), genes))
                        tmp_count = [[] for i in range(1 + len(filter_combs))]
                    t1 = time.time() - t0
                    print >> sys.stderr, '%i (last 100000 took %.2f secs)' % (counter, t1) 
                    t0 = time.time()

            counter += 1

            if options.is_bam:

                if line.is_unmapped:
                    continue

                if options.best_only and line.is_secondary:
                    continue

                #chrm = infile.getrname(line.tid).replace('chr', '')
                chrm = infile.getrname(line.tid)
                if chrm == 'chrM_rCRS':
                    chrm = 'chrM'
                pos = line.pos - 1
                broken = False

                #read_pos = line.positions --> alternative to code below
                read_pos = []
                for o in line.cigar:
                    if o[0] in [0, 2]:
                        read_pos.extend(range(pos, pos + o[1]))

                    if not o[0] in [1, 5]:
                        pos += o[1]

                try:
                    g = sp.unique(anno[chrm][read_pos])    
                except IndexError:
                    try:
                        read_pos = read_pos[(read_pos >= 0) & (read_pos < anno[chrm].shape[0])]
                        g = sp.unique(anno[chrm][read_pos])    
                    except:
                        continue

                g = g[g > 1]
                if g.shape[0] == 0:
                    continue

                ### resolve overlapping genes if we haven't masked them
                if not options.mask_gene_overlap and g.shape[0] > 1:
                    g = compress_g(g, idx2gene)
                tmp_count[0].extend(g)

                ### get validity for each filter
                if len(filter_names) > 0:
                    is_valid = sp.ones((len(filter_names), ), dtype = 'bool')
                    for idx, fn in enumerate(filter_names):
                        try:
                            is_valid[idx] = valid_after_filter(filters[idx][chrm], filter_types[idx], read_pos)
                        except KeyError:
                            continue
                    ### generate filter combination counts
                    for idx, comb in enumerate(filter_combs):
                        if sp.all(is_valid[comb]):
                            tmp_count[idx + 1].extend(g)
                    
            else:
                sl = line.strip().split('\t')
                if len(sl) < 9:
                    print >> sys.stderr, "ERROR: invalid SAM line\n%s" % line
                    sys.exit(1)

                (size, op) = (re.split('[^0-9]', sl[5])[:-1], re.split('[0-9]*', sl[5])[1:])
                size = [int(i) for i in size]

                #chrm = sl[2].replace('chr', '')
                chrm = sl[2]
                pos = int(sl[3]) - 1
                broken = False

                ## is unmapped ?
                if (int(sl[1]) & 4) == 4:
                    continue

                ## is secondary ?
                if options.best_only and (int(sl[1]) & 256 == 256):
                    continue

                for o in range(len(op)):
                    if op[o] in ['M', 'D']:
                        for p in range(size[o]):
                            try:
                                g = anno[chrm][pos + p]
                                if g > 1:
                                    tmp_count.append(g)
                                    break
                            except KeyError:
                                continue
                            except IndexError:
                                if chrm in ['chrM', 'M', 'chrM_rCRS']:
                                    continue
                                else:
                                    print >> sys.stderr, 'ERROR: %i exceeds length of %s' % (pos + p, chrm)
                    if broken:
                        break
                    if not op[o] in ['H', 'I']:
                        pos += size[o]
                           
        ### close file stream
        if not file == '-':
            infile.close()

        ### compress remaining counts
        for idx in range(len(tmp_count)):
            compressed_counts[idx].extend(compress_counts(sorted(tmp_count[idx]), genes))
        tmp_count = [[] for i in range(1 + len(filter_combs))]

        ### condense count lists
        print >> sys.stderr, "Sorting and condensing compressed list ..."
        compressed_counts = condense_compressed_counts(compressed_counts)

        ### resolve gene combinations
        for idx in range(len(compressed_counts)):
            extend_list = []
            for a in range(len(compressed_counts[idx]) -1, -1, -1):
                if len(idx2gene[compressed_counts[idx][a][0]]) > 1:
                    for g in idx2gene[compressed_counts[idx][a][0]]:
                        extend_list.append([gene2idx[tuple([g])], compressed_counts[idx][a][1]])
                    del compressed_counts[idx][a]
            compressed_counts[idx].extend(extend_list)
        compressed_counts = condense_compressed_counts(compressed_counts)

        ### remove gene IDs that encode combinations
        genes = [genes[i] for i in range(len(genes)) if len(idx2gene[genes[i]]) < 2]

        ### report counts to outfile
        if options.verbose:
            print >> sys.stderr, "Summarizing gene counts ..."

        for idx in range(len(compressed_counts)):
            if idx > 0 and (idx - 1) < len(filter_combs):
                comb_tag = '_'.join([filter_names[i] for i in filter_combs[idx - 1]])
            else:
                comb_tag = ''

            outfile = open(options.outfile + comb_tag, 'w')
            a = 0
            g = 0
            ### seek to first position that mapped to gene (0 means not gene found)
            while g < len(genes):
                while g < len(genes) and (a == len(compressed_counts[idx]) or genes[g] < compressed_counts[idx][a][0]):
                    print >> outfile, '%s\t0' % idx2gene[genes[g]][0]
                    if options.verbose and g % 100 == 0:
                        print >> sys.stderr, "%.2f / 100 percent \r" % (float(g) / len(genes) * 100),
                    g += 1
                while a < len(compressed_counts[idx]) and g < len(genes) and genes[g] == compressed_counts[idx][a][0]:
                    print >> outfile, '%s\t%i' % (idx2gene[genes[g]][0], compressed_counts[idx][a][1])
                    a += 1
                    g += 1
                    if options.verbose and g % 100 == 0:
                        print >> sys.stderr, "%.2f / 100 percent \r" % (float(g) / len(genes) * 100),
            outfile.close()

        if options.verbose:
            t1 = time.time() - time_total
            print >> sys.stderr, "\n... done - total run took %i secs." % t1

if __name__ == "__main__":
    main()
