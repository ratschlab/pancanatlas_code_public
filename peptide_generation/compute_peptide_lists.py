from __future__ import print_function

# * Python libraries
import sys
import cPickle
import os
import timeit
import ipdb
import argparse

# External libraries
import gff3
import h5py
import numpy as np
import Bio.SeqIO as BioIO

# Immunopeptide module

from immunopeptide.util_print import print_memory_diags, print_dict, \
    print_gene, print_segment_graph, print_splice_graph, print_line_feature, \
    print_ann_raw_file_lines, print_ann_file_lines, print_transcripts, print_reading_frames, \
    print_som_info, print_cds_begin_dict, print_annotation_result, print_paths
    
from immunopeptide.util_genomic import leq_strand, to_adj_list, to_adj_succ_list, \
    print_adj_list, translate_dna_to_peptide, has_stop_codon_initial, cross_peptide_result, \
    has_stop_codon_cross, complementary_seq, encode_chromosome

from immunopeptide.util_annotations import contained_ann_lines,  extract_parent_transcript, extract_attr, \
    filter_cds_lines, find_overlapping_cds, find_overlapping_cds_simple, preprocess_ann

from immunopeptide.util_misc import header_labels

from immunopeptide.util_io import save_list_csv, generate_output_fp, donor_list

from immunopeptide.util_splicegraph import genes_preprocess

from immunopeptide.util_mutations import apply_mutations, parse_somatic_mutation_info, parse_germline_mutation_info

from immunopeptide.util_expression import parse_gene_metadata_info, parse_gene_coverage, search_metadata_segmentgraph




## --------------------------- PEPTIDE GENERATION FUNCTIONS --------------------------------------------------------------------------------


# Optimized annotation code that does not loop over the annotation but uses the lookup structure that was built from the only initial pass
# over the GFF annotation file
def annotate_gene_opt(gene, seg_lookup_table, strain_idx_table, segment_expr_info,
                      cds_lookup_table, ref_seq, mut_seq, fa_ptr, donor_id, gene_norm_coverage, mutation_mode, som_pos_to_desc):
    peptide_type, mutation_type = header_labels(donor_id, mutation_mode)
    gene_name = gene.name
    cds_gene = gene.cds_gene
    nvertices = gene.nvertices
    vertex_succ_list = gene.vertex_succ_list
    vertex_order = gene.vertex_order
    gene.splicegraph.reading_frames = {}
    gene.splicegraph.peptides = {}

    # Renew peptides and reading frames associated with vertex wrt. donor ID
    for idx in range(nvertices):
        gene.splicegraph.reading_frames[idx] = set()
        gene.splicegraph.peptides[idx] = set()

    # Computation proceeds in vertex order in direction of read strand
    for idx in vertex_order:
        v_start = gene.splicegraph.vertices[0][idx]+1
        v_stop = gene.splicegraph.vertices[1][idx]
        cds_begins = gene.cds_begin_dict[idx]
        vertex_length = gene.vertex_len_dict[idx]

        # Skip de-generate exons that contain less than one codon
        if vertex_length < 3:
            continue

        # Initialize reading regions from the CDS transcript annotations
        for cds_begin in cds_begins:
            line_elems = cds_begin[1]["line_raw"].split()
            cds_strand = line_elems[6]
            assert(cds_strand == gene.strand)
            cds_phase = int(line_elems[7])

            if gene.strand == "-":
                cds_start = int(line_elems[4])
                n_trailing_bases = cds_start-cds_phase-v_start+1
                read_start_coord = cds_start-cds_phase
                read_stop_coord = v_start
            else:  # gene.strand=="+":
                cds_start = int(line_elems[3])
                n_trailing_bases = v_stop-cds_start-cds_phase+1
                read_start_coord = cds_start+cds_phase
                read_stop_coord = v_stop

            # The seeding is ill-formed if there are less than 3 bases on the vertex to be read
            if gene.strand == "+" and read_stop_coord-read_start_coord+1 < 3 or \
               gene.strand == "-" and read_start_coord-read_stop_coord+1 < 3:
                continue

            stop_codon_exist = has_stop_codon_initial(mut_seq, read_start_coord, read_stop_coord, gene.strand)

            # A read_phase of a vertex is defined as the number of overhanging bases with respect to a read frame sequence
            # A read frame is only propagated to the next vertex if a stop codon does not stop the read inside the vertex.
            if not stop_codon_exist:
                assert(n_trailing_bases >= 0)
                read_phase = n_trailing_bases % 3
                gene.splicegraph.reading_frames[idx].add(read_phase)

        vertex_seq_mut = mut_seq[gene.splicegraph.vertices[0, idx]:gene.splicegraph.vertices[1, idx]]
        vertex_seq_ref = ref_seq[gene.splicegraph.vertices[0, idx]:gene.splicegraph.vertices[1, idx]]
        n_read_frames = len(gene.splicegraph.reading_frames[idx])

        # Update the propagated read phases of the successor vertices.
        for prop_vertex in vertex_succ_list[idx]:

            prop_vertex_length = gene.splicegraph.vertices[1, prop_vertex]-gene.splicegraph.vertices[0, prop_vertex]
            prop_vertex_seq_mut = mut_seq[gene.splicegraph.vertices[0, prop_vertex]:gene.splicegraph.vertices[1, prop_vertex]]
            prop_vertex_seq_ref = ref_seq[gene.splicegraph.vertices[0, prop_vertex]:gene.splicegraph.vertices[1, prop_vertex]]

            # Skip de-generate exons that contain less than one codon
            if prop_vertex_length < 3:
                continue

            for read_frame in gene.splicegraph.reading_frames[idx]:
                seq_key = gene.name.strip()+"_"+str(idx)+"_"+str(prop_vertex)+"_"+str(read_frame)
                cross_peptide_mut, cross_peptide_ref, peptide_start_coord_v1, peptide_stop_coord_v1, \
                    peptide_start_coord_v2, peptide_stop_coord_v2, has_stop_codon = \
                    cross_peptide_result(read_frame, gene.strand,
                                         vertex_seq_ref, prop_vertex_seq_ref,
                                         vertex_seq_mut, prop_vertex_seq_mut,
                                         gene.splicegraph.vertices[:, idx],
                                         gene.splicegraph.vertices[:, prop_vertex])
                peptide_weight = 1.0/n_read_frames

                # If cross junction peptide has a stop-codon in it, we will not output it for simplicity, also the frame
                # will not be propagated because the read is truncated before it reaches the end of the exon.
                if not has_stop_codon:

                    # Write the variant gene into the FASTA FP together with the donor ID

                    if configs["gtex_normals"] or cross_peptide_mut!=cross_peptide_ref or donor_id=="ref":
                        peptide_type = "REFERENCE" if donor_id == "ref" else donor_id
                        header_line = " ".join([peptide_type, gene.name, gene.chr, gene.strand, "{:.3f}".format(gene_norm_coverage),
                                                "CROSSJUNCTION", mutation_type, "{:.3f}".format(peptide_weight)])
                        write_protein=True
                    else:
                        write_protein=False

                    if write_protein:
                        snp_descs = []
                        header_line += (" "+str(peptide_start_coord_v1)+" "+str(peptide_stop_coord_v1))
                        header_line += (" "+str(peptide_start_coord_v2)+" "+str(peptide_stop_coord_v2))

                        if mutation_mode is not None and (mutation_mode == "both" or mutation_mode == "somatic_only"):

                            assert(peptide_start_coord_v1 <= peptide_stop_coord_v1)
                            for gen_coord in range(peptide_start_coord_v1, peptide_stop_coord_v1+1):
                                if (donor_id, gen_coord) in som_pos_to_desc:
                                    snp_descs.append((gen_coord, som_pos_to_desc[(donor_id, gen_coord)]))

                            assert(peptide_start_coord_v2 <= peptide_stop_coord_v2)
                            for gen_coord in range(peptide_start_coord_v2, peptide_stop_coord_v2+1):
                                if (donor_id, gen_coord) in som_pos_to_desc:
                                    snp_descs.append((gen_coord, som_pos_to_desc[(donor_id, gen_coord)]))

                            header_line += " "+str(len(snp_descs))+" "

                            for gen_coord, snp_desc in snp_descs:

                                if not configs["gtex_normals"]:
                                    expr_level=search_metadata_segmentgraph(gene, gen_coord, seg_lookup_table, strain_idx_table, segment_expr_info, donor_id) 
                                else:
                                    expr_level=-1
                                    
                                header_line += " ".join([str(gen_coord), snp_desc["t_depth"], snp_desc["t_ref_count"], snp_desc["t_alt_count"], snp_desc["n_depth"],
                                                         snp_desc["n_ref_count"], snp_desc["n_alt_count"], str(expr_level)])

                        fa_ptr.write(">"+header_line+"\n")
                        peptide_str_pretty = "\n".join([cross_peptide_mut[i:i+100] for i in range(0, len(cross_peptide_mut), 100)])
                        fa_ptr.write(peptide_str_pretty+"\n\n")

                    n_leading_bases = (3-read_frame) % 3
                    propagated_frame = (prop_vertex_length-n_leading_bases) % 3
                    gene.splicegraph.reading_frames[prop_vertex].add(propagated_frame)

    gene.processed = True



# Generates a list of all peptides that are generated given the information in the GFF annotation file.
# At the same time it writes the generated peptides to the output FASTA file
def find_background_peptides(gene, seg_lookup_table, strain_idx_table, segment_expr_info, ref_seq, mut_seq, gene_ts_table, transcript_cds_table, output_fasta_fp, donor_id,
                             gene_norm_coverage, mutation_mode, som_pos_to_desc):
    gene_transcripts = gene_ts_table[gene.name]
    peptide_list = []
    peptide_type, mutation_type = header_labels(donor_id, mutation_mode)
    segment_graph = gene.segmentgraph

    # Generate a background peptide for every variant transcript
    for ts in gene_transcripts:

        # No CDS entries for transcript in annotation file...
        if ts not in transcript_cds_table:
            print("WARNING: Transcript not in CDS table")
            continue

        cds_list = transcript_cds_table[ts]

        # Reverse CDS in reverse strand transcription...
        if gene.strand.strip() == "-":
            cds_list = cds_list[::-1]

        cds_string_mutated = ""
        cds_string_unmutated = ""
        first_cds = True

        # Append transcribed CDS regions to the output
        for coord_left, coord_right, frameshift in cds_list:

            # Apply initial frameshift on the first CDS of the transcript
            if first_cds:
                if gene.strand.strip() == "+":
                    coord_left += frameshift
                else:
                    coord_right -= frameshift

            unmutated_nuc_seq = ref_seq[coord_left-1:coord_right]
            mutated_nuc_seq = mut_seq[coord_left-1:coord_right]

            # Accumulate new DNA sequence...
            if gene.strand.strip() == "+":
                cds_string_mutated += mutated_nuc_seq
                cds_string_unmutated += unmutated_nuc_seq
            elif gene.strand.strip() == "-":
                cds_string_mutated += complementary_seq(mutated_nuc_seq[::-1])
                cds_string_unmutated += complementary_seq(unmutated_nuc_seq[::-1])
            else:
                print("ERROR: Invalid strand...")
                sys.exit(1)

            first_cds = False

        # We will only write out variant transcripts that are differing from the
        # reference transcript, or in reference mode we always output.
        if cds_string_mutated != cds_string_unmutated or donor_id == "ref":
            aa_str_mutated = translate_dna_to_peptide(cds_string_mutated)
            peptide_list.append((aa_str_mutated, cds_list, ts, peptide_type, mutation_type))

    # Write a FASTA sequence to the output file for each variant peptide
    for peptide_str, coord_desc, transcript_id, sample_id, mutation_id in peptide_list:
        snp_descs = []
        n_gff_regions = len(coord_desc)
        header_line = " ".join([sample_id, gene.name, gene.chr, gene.strand, "{:.3f}".format(gene_norm_coverage), "BACKGROUND", mutation_id,
                                transcript_id, str(n_gff_regions)])
        for coding_region in coord_desc:
            region_start = coding_region[0]
            region_end = coding_region[1]

            if mutation_mode == "both" or mutation_mode == "somatic_only":
                for gen_coord in range(region_start, region_end+1):
                    if (donor_id, gen_coord) in som_pos_to_desc:
                        snp_descs.append((gen_coord, som_pos_to_desc[(donor_id, gen_coord)]))

            header_line += (" "+str(region_start)+" "+str(region_end))

        if mutation_mode is not None and (mutation_mode == "both" or mutation_mode == "somatic_only"):
            header_line += " "+str(len(snp_descs))+" "

            for gen_coord, snp_desc in snp_descs:
                
                if not configs["gtex_normals"]:
                    expr_level=search_metadata_segmentgraph(gene, gen_coord, seg_lookup_table, strain_idx_table, segment_expr_info, donor_id)
                else:
                    expr_level=-1

                header_line += " ".join([str(gen_coord), snp_desc["t_depth"], snp_desc["t_ref_count"], snp_desc["t_alt_count"], snp_desc["n_depth"],
                                         snp_desc["n_ref_count"], snp_desc["n_alt_count"], str(expr_level)])

        output_fasta_fp.write(">"+header_line+"\n")
        peptide_str_pretty = "\n".join([peptide_str[i:i+100] for i in range(0, len(peptide_str), 100)])
        output_fasta_fp.write(peptide_str_pretty+"\n\n")
    gene.processed = True


# ----------------------------------------- PEPTIDE LIST I/O FUNCTIONS -----------------------------------------------------------------------------------


# Writes out personalized peptides to file, also accepts special argument 'reference', which performs the same
#  analysis but without applying mutations, but only against the reference genome.
def write_peptides_donor(seq_dict, som_info, som_pos_to_desc, germ_info, gexpr_info, graph_data, seg_lookup_table, strain_idx_table, segment_expr_info,
                         cds_lookup_table, gene_to_transcript_table, transcript_to_cds_table, output_fasta_fp, donor_id="ref", mutation_mode=None):
    assert(not (donor_id == "ref" and mutation_mode is not None))

    # In an initial pseudo-run, the peptides emitted from the background are written out...
    if donor_id == "ref":
        print("Generating reference sequences...")
        donor_seq_dict = seq_dict
    else:
        start_time = timeit.default_timer()

        if donor_id not in som_info:
            print("WARNING: {} has no somatic mutation records, skipping".format(donor_id))
            return

        if mutation_mode == "both":
            donor_seq_dict = apply_mutations(seq_dict, apply_somatic=True, som_info=som_info[donor_id],
                                             apply_germline=True, germ_info=germ_info[donor_id])
        elif mutation_mode == "somatic_only":
            donor_seq_dict = apply_mutations(seq_dict, apply_somatic=True, som_info=som_info[donor_id],
                                             apply_germline=False, germ_info=germ_info[donor_id])
        elif mutation_mode == "germline_only":
            donor_seq_dict = apply_mutations(seq_dict, apply_somatic=False, som_info=som_info[donor_id],
                                             apply_germline=True, germ_info=germ_info[donor_id])
        else:
            print("ERROR: Invalid mutation mode...")
            sys.exit(1)

        end_time = timeit.default_timer()
        print("MUTATION APPLY TIME: {:.3f} seconds".format(end_time-start_time))

    cum_time_bg = 0.0
    cum_time_cross = 0.0
    n_genes_processed = 0

    # Process each of the interesting genes independently...
    for gene_idx in np.arange(graph_data.shape[0]):

        gene = graph_data[gene_idx]
        gene_norm_coverage = np.nan if donor_id == "ref" else gexpr_info[(gene.name.strip(), donor_id.strip())]

        # Genes not contained in the annotation...
        if gene.name not in cds_lookup_table or gene.name not in gene_to_transcript_table:
            gene.processed = False
            continue

        # Logging...
        if gene_idx % 1000 == 0:
            print("Processing gene {}/{}: {}".format(gene_idx, graph_data.shape[0], gene.name))
            print("BG TIME PER GENE: {:.3f}".format(cum_time_bg/(gene_idx+1)))
            print("CROSS TIME PER GENE: {:.3f}".format(cum_time_cross/(gene_idx+1)))
            print("Number of completed genes: {}".format(n_genes_processed))

        gene_length = gene.stop-gene.start
        ref_seq = seq_dict[gene.chr.strip()]
        mut_seq = donor_seq_dict[gene.chr.strip()]
        start_time = timeit.default_timer()

        if not configs["gtex_normals"]:
            find_background_peptides(gene,seg_lookup_table, strain_idx_table, segment_expr_info, ref_seq,mut_seq,gene_to_transcript_table,transcript_to_cds_table,
                                     output_fasta_fp,donor_id,gene_norm_coverage, mutation_mode, som_pos_to_desc)

        end_time = timeit.default_timer()
        cum_time_bg += end_time-start_time
        start_time = timeit.default_timer()
        annotate_gene_opt(gene, seg_lookup_table, strain_idx_table, segment_expr_info, cds_lookup_table, ref_seq, mut_seq, output_fasta_fp, donor_id,
                          gene_norm_coverage, mutation_mode, som_pos_to_desc)
        end_time = timeit.default_timer()
        cum_time_cross += end_time-start_time
        n_genes_processed += 1


# Write the resulting peptides in FASTA format
def write_peptides(som_info, som_pos_to_desc, germ_info, gexpr_info, donor_list, input_instance, testmode=False, parallel_rank=None, configs=None):
    seq_file_path=configs["sequence_path"]
    graph_path=configs["graph_path"]
    ann_path=configs["annotation_path"]
    graph_metadata_file=configs["graph_metadata_path"]
    graph_fp=open(graph_path, 'r')
    seq_dict = {}
    start_time = timeit.default_timer()
    interesting_chr = map(str, range(1, 23))+["X", "Y", "MT"]

    # Read reference genome for standard chromosomes
    for record in BioIO.parse(seq_file_path, "fasta"):
        if record.id in interesting_chr:
            seq_dict[record.id] = str(record.seq).strip()

    end_time = timeit.default_timer()
    print("SEQUENCE PARSE TIME: {:.3f} seconds".format(end_time-start_time))
    print_memory_diags()
    start_time = timeit.default_timer()
    graph_obj = cPickle.load(graph_fp)
    graph_fp.close()
    graph_data = graph_obj[0]
    graph_meta = graph_obj[1]
    end_time = timeit.default_timer()
    print("SPLICE GRAPH LOAD TIME: {:.3f} seconds".format(end_time-start_time))
    print_memory_diags()
    start_time = timeit.default_timer()

    h5f = h5py.File(graph_metadata_file, 'r')
    print("LOADED HDF5 DATASET...")
    print_memory_diags()
    seg_lookup_table, strain_idx_table, segment_expr_info = parse_gene_metadata_info(h5f, donor_list)
    end_time = timeit.default_timer()
    print("SEGMENT GRAPH META-DATA LOAD TIME: {:.3f} seconds".format(end_time-start_time))
    print_memory_diags()
    start_time = timeit.default_timer()
    ann_obj = gff3.Gff3()
    ann_obj.parse(ann_path)
    end_time = timeit.default_timer()
    print("GFF STRUCTURE PARSE TIME: {:.3f} seconds".format(end_time-start_time))
    print_memory_diags()
    line_view = ann_obj.lines
    start_time = timeit.default_timer()
    cds_lookup_table, gene_to_transcript_table, transcript_to_cds_table = preprocess_ann(line_view)
    end_time = timeit.default_timer()
    print("BUILD LOOKUP STRUCTURE TIME: {:.3f} seconds".format(end_time-start_time))
    print_memory_diags()
    cum_time = 0.0
    cum_time_bg = 0.0
    cum_time_cross = 0.0
    n_gene_unknown_cds = 0
    n_gene_unknown_transcript = 0
    start_time = timeit.default_timer()
    genes_preprocess(graph_data, cds_lookup_table)
    end_time = timeit.default_timer()
    print("GENES PREPROCESSING TIME: {:.3f} seconds".format(end_time-start_time))
    assert(len(donor_list) == 1)

    # The first job gets to write out the reference sequence into its file.
    if parallel_rank == 0:
        output_fasta_fp = generate_output_fp(input_instance, "ref", "background", configs)
        write_peptides_donor(seq_dict, som_info, som_pos_to_desc, germ_info, gexpr_info, graph_data, seg_lookup_table, strain_idx_table, segment_expr_info,
                             cds_lookup_table, gene_to_transcript_table, transcript_to_cds_table, output_fasta_fp, donor_id="ref")
        output_fasta_fp.close()

    # Personalized analysis
    for donor_idx, donor_id in enumerate(donor_list):
        print("Processing donor: {} {}/{}".format(donor_id, donor_idx+1, len(donor_list)))
        start_time = timeit.default_timer()
        output_fasta_fp = generate_output_fp(input_instance, donor_id, "somatic_only", configs)
        write_peptides_donor(seq_dict, som_info, som_pos_to_desc, germ_info, gexpr_info, graph_data, seg_lookup_table, strain_idx_table, segment_expr_info,
                             cds_lookup_table, gene_to_transcript_table, transcript_to_cds_table, output_fasta_fp,
                             donor_id=donor_id, mutation_mode="somatic_only")
        output_fasta_fp.close()
        output_fasta_fp = generate_output_fp(input_instance, donor_id, "germline_only", configs)
        write_peptides_donor(seq_dict, som_info, som_pos_to_desc, germ_info, gexpr_info, graph_data, seg_lookup_table, strain_idx_table, segment_expr_info,
                             cds_lookup_table, gene_to_transcript_table, transcript_to_cds_table, output_fasta_fp,
                             donor_id=donor_id, mutation_mode="germline_only")
        output_fasta_fp.close()
        output_fasta_fp = generate_output_fp(input_instance, donor_id, "both", configs)
        write_peptides_donor(seq_dict, som_info, som_pos_to_desc, germ_info, gexpr_info, graph_data, seg_lookup_table, strain_idx_table, segment_expr_info,
                             cds_lookup_table, gene_to_transcript_table, transcript_to_cds_table, output_fasta_fp,
                             donor_id=donor_id, mutation_mode="both")
        output_fasta_fp.close()
        end_time = timeit.default_timer()
        cum_time += (end_time-start_time)
        print("AVERAGE PROCESSING TIME PER DONOR: {:.3f} seconds".format(cum_time/(donor_idx+1)))

    h5f.close()


#  ------------------------------------------ SCRIPT ENTRY POINT ------------------------------------------------------------------

def compute_peptide_lists(configs):
    ''' Main entry point for the script'''
    sys.path.append(configs["spladder_path"])
    parallel_rank = configs["parallel_rank"]
    print("Parallel rank: {}".format(parallel_rank))
    log_path = configs["log_path"]

    if configs["batch_mode"]:
        sys.stdout = open(os.path.join(log_path, "{}_stdout.log".format(parallel_rank)), 'w', 0)
        sys.stderr = open(os.path.join(log_path, "{}_stderr.log".format(parallel_rank)), 'w', 0)

    input_instance = configs["run_id"]
    don_list_full = donor_list(configs)
    assert(len(don_list_full) == 63)
    don_list = [don_list_full[parallel_rank]]
    start_time = timeit.default_timer()
    germ_info = parse_germline_mutation_info(don_list, configs)
    end_time = timeit.default_timer()
    proc_time = end_time-start_time
    print("Parse germline mutation time: {:.2f}".format(proc_time))
    print_memory_diags()
    start_time = timeit.default_timer()
    som_info, som_pos_to_desc = parse_somatic_mutation_info(don_list, configs)
    end_time = timeit.default_timer()
    proc_time = end_time-start_time
    print("Parse somatic mutation time: {:.2f}".format(proc_time))
    print_memory_diags()
    start_time = timeit.default_timer()
    gexpr_info = parse_gene_coverage(don_list, configs)
    end_time = timeit.default_timer()
    proc_time = end_time-start_time
    print("Parse gene coverage time: {:.2f}".format(proc_time))
    print_memory_diags()
    write_peptides(som_info, som_pos_to_desc, germ_info, gexpr_info, don_list, input_instance, parallel_rank=parallel_rank, configs=configs)


if __name__ == "__main__":

    parser=argparse.ArgumentParser()
    
    # ----- Substitute the local code/data base path here-----
    
    INPUT_DATA_PATH=None
    OUTPUT_DATA_PATH=None
    CODE_PATH=None

    # -------------------------------------------------------------------------------------------------------------------------------------

    # Input / library paths
    
    ANNOTATION_PATH=os.path.join(INPUT_DATA_PATH,"annotation","gencode.v19.annotation.hs37d5_chr.gff")
    
    # One of {genes_graph_conf3.merge_graphs.validated.pickle, genes_graph_gtex_conf2.merge_graphs.pickle}, the first
    # one should be used for generating TCGA peptides, the second one for generating GTEX background peptides

    GTEX_GRAPH_PATH=os.path.join(INPUT_DATA_PATH,"graph","genes_graph_gtex_conf2.merge_graphs.pickle")
    TCGA_GRAPH_PATH=os.path.join(INPUT_DATA_PATH,"graph","genes_graph_conf3.merge_graphs.validated.pickle")
    GRAPH_PATH=TCGA_GRAPH_PATH

    # One of {genes_graph_conf3.merge_graphs.validated.count.h5, genes_graph_gtex_conf2.merge_graphs.validated.count.h5}, the first
    # one should be used for generating TCGA peptides, the second one for generating GTEX background peptides
    
    GTEX_GRAPH_METADATA_PATH=os.path.join(INPUT_DATA_PATH,"graph_metadata","genes_graph_gtex_conf2.merge_graphs.validated.count.h5")
    TCGA_GRAPH_METADATA_PATH=os.path.join(INPUT_DATA_PATH,"graph_metadata","genes_graph_conf3.merge_graphs.validated.count.h5")
    GRAPH_METADATA_PATH=TCGA_GRAPH_METADATA_PATH

    # Sequence
    SEQ_PATH=os.path.join(INPUT_DATA_PATH,"sequence","genome.fa")

    # Somatic mutations
    SOM_PATH=os.path.join(INPUT_DATA_PATH,"somatic_variants","pancan.merged.v0.2.6.PUBLIC.maf")

    # Germline mutations
    GERMLINE_PATH=os.path.join(INPUT_DATA_PATH,"germline_variants","mergedfiles_clean_stringentfilter.h5")

    # Gene expression information
    GEXPR_PATH=os.path.join(INPUT_DATA_PATH,"gene_expression","expression_for_immunology.h5")

    # List of donor sample IDs to be generated
    DONOR_PATH=os.path.join(INPUT_DATA_PATH,"sample_lists","donor_cptac_rerun_181105.tsv")

    # Spladder dependency
    SPLADDER_PATH=os.path.join(CODE_PATH,"spladder","python")

    # ---------------------------------------------------------------------------------------------------------------------------------------
    
    # Output paths
    
    OUTPUT_PATH=os.path.join(OUTPUT_DATA_PATH,"peptide_output")
    LOG_PATH=os.path.join(OUTPUT_DATA_PATH,"peptide_logs")

    # Paths
    parser.add_argument("--annotation_path", default=ANNOTATION_PATH, help="File storing the canonical transcript annotations")
    parser.add_argument("--graph_path", default=GRAPH_PATH, help="File storing the splicing graph to be processed")
    parser.add_argument("--graph_metadata_path", default=GRAPH_METADATA_PATH, help="File storing the graph meta-data")
    parser.add_argument("--sequence_path", default=SEQ_PATH, help="File storing the reference genome sequence")
    parser.add_argument("--somatic_mutation_path", default=SOM_PATH, help="File storing the somatic variants to be attached")
    parser.add_argument("--germline_mutation_path", default=GERMLINE_PATH, help="File storing the germline variants to be attached")
    parser.add_argument("--gene_expr_path", default=GEXPR_PATH, help="File storing the gene expressions")
    parser.add_argument("--donor_list_path", default=DONOR_PATH, help="File specifying the donor sample IDs to be processed")
    parser.add_argument("--spladder_path", default=SPLADDER_PATH, help="Location of spladder dependency")
    parser.add_argument("--output_path", default=OUTPUT_PATH, help="Output base paths where peptide lists are stored")
    parser.add_argument("--log_path", default=LOG_PATH, help="Logging location")

    # Arguments
    parser.add_argument("--gtex_normals", default=False, action="store_true", help="GTEX normal mode?")
    parser.add_argument("--run_id", default="peptidegen_run", help="Description of the run, output files will be labeled with it")
    parser.add_argument("--parallel_rank", default=0, type=int, help="Parallel rank that is used in batch processing")
    parser.add_argument("--batch_mode", default=False, action="store_true", help="Batch mode?")

    args=parser.parse_args()
    configs=vars(args)

    compute_peptide_lists(configs)

