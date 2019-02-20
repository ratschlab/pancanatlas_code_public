''' Utility functions concerned with parsing of annotation file'''

import sys

from immunopeptide.util_genomic import leq_strand

# Returns the lines of an annotation file that are strictly inclusive of gene coordinate
#  bounds
def contained_ann_lines(start_coord, stop_coord, gene_name, ann_lines):
    ret_list = []
    for idx, line in enumerate(ann_lines):
        if line["line_type"].strip() in ["directive", "comment"]:
            continue
        raw_line = line["line_raw"]
        line_elems = raw_line.split()
        elem_start = int(line_elems[3])
        elem_end = int(line_elems[4])
        if elem_start >= start_coord and elem_end <= stop_coord:
            ret_list.append(line)
    return ret_list

# Returns the parent transcript from a line
def extract_parent_transcript(cds_line):
    attr_elems = cds_line.split()[8].split(";")
    parent_attr = filter(lambda attr_elem: "Parent=" in attr_elem, attr_elems)
    assert(len(parent_attr) == 1)
    parent_ts = parent_attr[0].split("=")[1].strip()
    return parent_ts

# Extract GeneID from line
def extract_attr(cds_line, attr_key):
    attr_elems = cds_line.split()[8].split(";")
    attr = filter(lambda attr_elem: attr_elem.startswith(attr_key+"="), attr_elems)
    assert(len(attr) == 1)
    attr_val = attr[0].split("=")[1].strip()
    return attr_val




# Returns the CDS lines that define a start CDS, i.e. the ones first in genomic
#  coordinates for all the CDS, only returns CDS that are associated with the
#  <target_gene>
def filter_cds_lines(ann_lines, strand_mode, target_gene):
    target_transcripts = set()
    for line in ann_lines:
        if line["line_type"].strip() in ["directive", "comment"]:
            continue
        raw_line = line["line_raw"]
        raw_line_elems = raw_line.split()
        feature_type = raw_line_elems[2]
        if "transcript" not in feature_type and "mRNA" not in feature_type:
            continue
        gene_id = extract_attr(raw_line, "geneID")
        if not gene_id == target_gene:
            continue
        transcript_id = extract_attr(raw_line, "ID")
        target_transcripts.add(transcript_id)
    print("Number of transcripts/mRNA: {}".format(len(target_transcripts)))
    smallest_cds = {}

    # Traverse annotation lines
    for line in ann_lines:
        if line["line_type"].strip() in ["directive", "comment"]:
            continue
        raw_line = line["line_raw"]
        raw_line_elems = raw_line.split()
        feature_type = raw_line_elems[2]
        if "CDS" not in feature_type:
            continue
        parent_ts = extract_parent_transcript(raw_line)

        if parent_ts not in target_transcripts:
            continue

        cds_start = int(raw_line_elems[3]) if strand_mode == "+" else int(raw_line_elems[4])

        if parent_ts not in smallest_cds or leq_strand(cds_start, smallest_cds[parent_ts][0], strand_mode):
            smallest_cds[parent_ts] = (cds_start, line)

    relevant_cds_lines = map(lambda pair: pair[1], smallest_cds.values())
    print("Number of relevant CDS begin: {}".format(len(relevant_cds_lines)))
    return relevant_cds_lines



# Find overlapping CDS regions
# cds_lines: Set of CDS lines in annotation file
# v_start: Start coordinate of exon
# v_stop: Stop coordinate of exon
# strand: +/- with + reading towards positive coordinates and - reading towards negative coordinates
def find_overlapping_cds(cds_lines, v_start, v_stop, strand):
    if strand == "+":
        read_begin_idx = 3
    elif strand == "-":
        read_begin_idx = 4
    else:
        print("ERROR: Read strand not given")
        sys.exit(1)
    return filter(lambda line: int(line["line_raw"].split()[read_begin_idx]) >= v_start
                  and int(line["line_raw"].split()[read_begin_idx]) <= v_stop, cds_lines)


# Find overlapping CDS given a list of CDS starts
def find_overlapping_cds_simple(v_start, v_stop, cds_begins):
    return filter(lambda cds_begin: cds_begin[0] >= v_start and cds_begin[0] <= v_stop, cds_begins)



# Pre-processed the annotation file and builds a lookup structure that can be used to retrieve
#  the CDS beginnings when looping over the genes
def preprocess_ann(ann_lines):
    transcript_to_gene_dict = {}
    gene_to_transcript_dict = {}
    transcript_to_cds_dict = {}
    transcript_cds_begin_dict = {}
    gene_cds_begin_dict = {}
    for idx, line in enumerate(ann_lines):
        if line["line_type"].strip() in ["directive", "comment"]:
            continue
        raw_line = line["line_raw"]
        raw_line_elems = raw_line.split()
        feature_type = raw_line_elems[2]
        if feature_type in ["transcript", "mRNA"]:
            gene_id = extract_attr(raw_line, "geneID")
            gene_type = extract_attr(raw_line, "gene_type")
            transcript_type = extract_attr(raw_line, "transcript_type")
            if gene_type.strip() != "protein_coding" or transcript_type.strip() != "protein_coding":
                continue
            transcript_id = extract_attr(raw_line, "ID")
            assert(transcript_id not in transcript_to_gene_dict)
            transcript_to_gene_dict[transcript_id] = gene_id
            if gene_id not in gene_to_transcript_dict:
                gene_to_transcript_dict[gene_id] = set()
            gene_to_transcript_dict[gene_id].add(transcript_id)
        elif feature_type == "CDS":
            parent_ts = extract_parent_transcript(raw_line)
            strand_mode = raw_line_elems[6]
            cds_left = int(raw_line_elems[3])
            cds_right = int(raw_line_elems[4])
            frameshift = int(raw_line_elems[7])
            if parent_ts not in transcript_to_cds_dict:
                transcript_to_cds_dict[parent_ts] = []
            transcript_to_cds_dict[parent_ts].append((cds_left, cds_right, frameshift))
            cds_start = cds_left if strand_mode == "+" else cds_right
            if parent_ts not in transcript_cds_begin_dict or leq_strand(cds_start, transcript_cds_begin_dict[parent_ts][0], strand_mode):
                transcript_cds_begin_dict[parent_ts] = (cds_start, line)
    n_transcripts = len(transcript_to_gene_dict.keys())
    for ts_key in transcript_to_gene_dict:
        target_gene = transcript_to_gene_dict[ts_key]
        if target_gene not in gene_cds_begin_dict:
            gene_cds_begin_dict[target_gene] = []
        if ts_key in transcript_cds_begin_dict:
            cds_start = transcript_cds_begin_dict[ts_key]
            gene_cds_begin_dict[target_gene].append(cds_start)
    for ts_key in transcript_to_cds_dict:
        cds_list = transcript_to_cds_dict[ts_key]
        transcript_to_cds_dict[ts_key] = sorted(cds_list, key=lambda coordpair: coordpair[0])
    return gene_cds_begin_dict, gene_to_transcript_dict, transcript_to_cds_dict
