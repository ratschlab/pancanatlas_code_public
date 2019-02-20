
''' Utility functions for printing internal data structures'''

import psutil
import os

from immunopeptide.util_genomic import to_adj_list

# Pretty prints a splicing graph
# splice_graph: <spladder.modules.classes.Splicegraph>
def print_splice_graph(splice_graph):
    nvertices = splice_graph.vertices.shape[1]
    print("Number of vertices: {}".format(nvertices))
    print("Number of edges: {}".format(splice_graph.edges.shape))
    print("Number of terminals: {}".format(splice_graph.terminals.shape))
    print("VERTEX ITERATION: ")
    for idx in range(nvertices):
        print("Vertex begin: {}".format(splice_graph.vertices[0, idx]))
        print("Vertex end: {}".format(splice_graph.vertices[1, idx]))
        print("Terminal map 0: {}".format(splice_graph.terminals[0, idx]))
        print("Terminal map 1: {}".format(splice_graph.terminals[1, idx]))
    print("EDGE ITERATION: ")
    adj_list = to_adj_list(splice_graph.edges)
    print(str(adj_list))


# Pretty prints a segment graph
# seg_graph: <spladder.modules.classes.Segmentgraph>
def print_segment_graph(seg_graph):
    pass

def print_memory_diags():
    ''' Print memory diagnostics including the active resident set size'''
    process = psutil.Process(os.getpid())
    print("Memory usage: {:.3f} GB".format(process.memory_info().rss/1000000000.0))

# Pretty-prints a dictionary
def print_dict(in_dict):
    for key in sorted(in_dict.keys()):
        print("{}: {}".format(key, in_dict[key]))

# Pretty prints information about a gene
# in_gene: <spladder.modules.classes.gene> describing a gene
def print_gene(in_gene, short=False):
    if short:
        print("Name: {}".format(in_gene.name))
    else:
        print("Name: {}".format(in_gene.name))
        print("Start: {}".format(in_gene.splicegraph.vertices.min()))
        print("Stop: {}".format(in_gene.splicegraph.vertices.max()))
        print("Gene length: {}".format(in_gene.stop-in_gene.start))
        print("Number of vertices: {}".format(in_gene.splicegraph.vertices.shape[1]))
        print("Read strand: {}".format(in_gene.strand))
        print("Number of transcripts: {}".format(len(in_gene.transcripts)))
        print("Number of exons: {}".format(len(in_gene.exons)))
        print("Gene type: {}".format(in_gene.gene_type))
        print("Gene source: {}".format(in_gene.source))
        print("Gene chromosome: {}".format(in_gene.chr))
        print("Has alt stand/end?: {}".format(in_gene.is_alt))
        print("Has alt splicing?: {}".format(in_gene.is_alt_spliced))



# Pretty prints a line feature
# in_line: A line feature of a GFF3 parsed structure
def print_line_feature(in_line):
    print(str(in_line["line_raw"]))

def print_ann_raw_file_lines(ann_file_lines):
    print("======== ANNOTATION FILE LINES =========")
    for idx, line in enumerate(ann_file_lines):
        print("LINE INDEX: {}".format(idx))
        print_line_feature(line)
    print("========= ANNOTATION FILE END ==========")

# Prints annotation file raw lines
def print_ann_file_lines(ann_file_lines):
    print("======== ANNOTATION FILE LINES =========")
    line_types = []
    print("Number of lines: {}".format(len(ann_file_lines)))
    for idx, line in enumerate(ann_file_lines):
        print("LINE INDEX: {}".format(idx))
        print_line_feature(line)
        line_elems = line["line_raw"].split()
        line_type = line_elems[2]
        elem_start = int(line_elems[3])
        elem_end = int(line_elems[4])
        assert(elem_start <= elem_end)
        line_types.append(line_type)
    line_types = set(line_types)
    print("Line types: {}".format(line_types))
    print("========= ANNOTATION FILE END ==========")

# Prints gene transcripts
def print_transcripts(gene):
    print("Gene transcripts BEGIN")
    for idx, tscript in enumerate(gene.transcripts):
        print("{} :\n{}".format(tscript, gene.exons[idx]))
    print("Gene transcripts END")


# Given a splice graph prints out the annotated reading frames of all vertices
def print_reading_frames(splicegraph, strand):
    nvertices = splicegraph.vertices.shape[1]
    sep_str = "->" if strand == "+" else "<-"
    for idx in range(nvertices):
        print("VERTEX {}: {}{}{}, READING FRAMES: {}".format(idx, splicegraph.vertices[0, idx], sep_str,
                                                             splicegraph.vertices[1, idx], splicegraph.reading_frames[idx]))


# Overview of somatic mutation information
def print_som_info(som_table):
    for k, v in som_table.iteritems():
        print("Somatic mutations for donor: {}".format(k))
        print("Number of mutations: {}".format(len(v)))


def print_cds_begin_dict(cds_begin_dict):
    for k in cds_begin_dict.keys():
        print("{}: ".format(k))
        cds_begins = cds_begin_dict[k]
        for idx, cds_begin in enumerate(cds_begins):
            print("{}: {}".format(idx, cds_begin[0]))



def print_annotation_result(gene):
    print("GENE ANNOTATION RESULT")
    if gene.processed:
        nvertices = gene.splicegraph.vertices.shape[1]
        for idx in range(nvertices):
            print("Vertex: {} [{},{}]".format(idx, gene.splicegraph.vertices[0, idx], gene.splicegraph.vertices[1, idx]))
            print("Reading frames: {}".format(gene.splicegraph.reading_frames[idx]))
            print("Number of generated peptides: {}".format(len(gene.splicegraph.peptides[idx])))
            for jdx, peptide in enumerate(gene.splicegraph.peptides[idx]):
                print("Peptide {}: {}\nweight: {}, v_from: {}, v_to: {}, start_coord: {}, stop_coord: {}".format(jdx, peptide[0],
                                                                                                                 peptide[1], idx,
                                                                                                                 peptide[2], peptide[3],
                                                                                                                 peptide[4]))
    else:
        print("Gene not processed, skipping...")


def print_paths(ann_path, graph_path, seq_path, som_path, donor_path):
    print("Annotation path: {}".format(ann_path))
    print("Graph path: {}".format(graph_path))
    print("Sequence path: {}".format(seq_path))
    print("Somatic mutation path: {}".format(som_path))
    print("Cancer type path: {}".format(donor_path))
