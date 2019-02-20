''' Utility functions related to splice-graph'''

import numpy as np

from immunopeptide.util_genomic import to_adj_succ_list
from immunopeptide.util_annotations import find_overlapping_cds_simple

# Pre-process gene structures to aid fast generation of cross-junction peptides
def genes_preprocess(genes, cds_lookup_table):

    for gene_idx in np.arange(genes.shape[0]):
        gene = genes[gene_idx]
        gene.from_sparse()
        assert(gene.strand in ["+", "-"])
        assert(len(gene.transcripts) == len(gene.exons))
        gene_name = gene.name

        # Ignore genes that are not annotated...
        if gene_name not in cds_lookup_table.keys():
            continue

        gene.cds_gene = cds_lookup_table[gene_name]
        gene.nvertices = gene.splicegraph.vertices.shape[1]
        gene.vertex_succ_list = to_adj_succ_list(gene.splicegraph.edges, gene.splicegraph.vertices, gene.strand)
        if gene.strand == "+":
            gene.vertex_order = np.argsort(gene.splicegraph.vertices[0, :])
        else:  # gene.strand=="-"
            gene.vertex_order = np.argsort(gene.splicegraph.vertices[1, :])[::-1]
        gene.cds_begin_dict = {}
        gene.vertex_len_dict = {}
        for idx in gene.vertex_order:
            v_start = gene.splicegraph.vertices[0][idx]+1
            v_stop = gene.splicegraph.vertices[1][idx]
            gene.cds_begin_dict[idx] = find_overlapping_cds_simple(v_start, v_stop, cds_lookup_table[gene_name])
            gene.vertex_len_dict[idx] = gene.splicegraph.vertices[1, idx]-gene.splicegraph.vertices[0, idx]
        gene.to_sparse()
