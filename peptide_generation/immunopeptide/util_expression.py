''' Utility functions concerned with expression of genomic regions'''

import h5py
import numpy as np

# Constructs look-up tables to lookup the correct column in matrix for a particular donor and
#  splits the global segments matrix into smaller chunks for each geneID to enable faster look-up
#  during peptide emission.
def parse_gene_metadata_info(h5f, donor_list):
    strain_expr_info = h5f["/strains"]
    segment_expr_info = h5f["/segments"]
    assert(strain_expr_info.size == segment_expr_info.shape[1])
    strain_idx_table = {}

    for strain_idx in np.arange(strain_expr_info.size):
        strain_id = strain_expr_info[strain_idx][0:12]

        if strain_id in donor_list:
            strain_idx_table[strain_id] = strain_idx

    gene_names = h5f["/gene_names"]
    gene_ids_segs = h5f["/gene_ids_segs"]
    assert(gene_ids_segs.size == segment_expr_info.shape[0])
    seg_lookup_table = {}

    for seg_idx in np.arange(gene_ids_segs.shape[0]):
        gene_id = gene_names[gene_ids_segs[seg_idx, 0]]

        if gene_id not in seg_lookup_table:
            seg_lookup_table[gene_id] = []

        seg_lookup_table[gene_id].append(seg_idx)

    return (seg_lookup_table, strain_idx_table, segment_expr_info)



# Parses the gene coverage information to propagate uncertainty into the peptide output
def parse_gene_coverage(donor_list, configs):
    gene_coverage_table = {}
    gexpr_path=configs["gene_expr_path"]

    with h5py.File(gexpr_path, 'r') as h5f:
        gene_ids = np.array(h5f["geneid"])
        donor_ids = np.array(h5f["donorid"])
        sel_lst = []
        run_idx_map = {}
        run_idx = 0

        for idx, elem in np.ndenumerate(donor_ids):
            tcga_key = '-'.join(donor_ids[idx[0]].split('-')[0:3])
            if tcga_key in donor_list:
                sel_lst.append(idx[0])
                run_idx_map[tcga_key] = run_idx
                run_idx += 1

        count_arr = np.array(h5f["counts"][:, sel_lst])

    for donor in donor_list:
        donor_col_idx = run_idx_map[donor]

        for row_idx in np.arange(count_arr.shape[0]):
            count = count_arr[row_idx, donor_col_idx]
            gene_coverage_table[(gene_ids[row_idx].strip(), donor)] = count

    return gene_coverage_table


# Given a segment graph and a 1-based genomic coordinate position of a somatic mutation, returns the expression
#  level of the one segment that contains the position. Performs a linear search in the segments contained in
#  the graph. Could be optimized to do a binary search, let's see if it is a bottleneck.
def search_metadata_segmentgraph(gene, gen_coord, seg_lookup_table, strain_idx_table, segment_expr_info, donor_id):
    gene_name = gene.name
    segment_graph = gene.segmentgraph
    seg_idxs = seg_lookup_table[gene_name]
    col_idx = strain_idx_table[donor_id]
    assert(len(seg_idxs) == segment_graph.segments.shape[1])

    # We use that in a segment graph, the vertices are not overlapping, so once we have found the segment,
    # we can look-up its expression value and return directly.
    for per_gene_idx, global_seg_idx in enumerate(seg_idxs):
        lbound = segment_graph.segments[0, per_gene_idx]+1
        rbound = segment_graph.segments[1, per_gene_idx]

        if gen_coord >= lbound and gen_coord <= rbound:
            gene_expr_entry = float(segment_expr_info[global_seg_idx, col_idx])
            return gene_expr_entry

    assert(False)
    return np.nan



