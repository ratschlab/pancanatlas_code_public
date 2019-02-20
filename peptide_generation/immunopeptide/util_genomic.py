
# Returns strand-sensitive order between two genomic coordinates
def leq_strand(coord1, coord2, strand_mode):
    if strand_mode == "+":
        return coord1 <= coord2
    else:  # strand_mode == "-"
        return coord1 >= coord2


# Converts a binary adjacency matrix to a list of directed edges
def to_adj_list(adj_matrix):
    adj_list = []
    assert(adj_matrix.shape[0] == adj_matrix.shape[1])
    for idx in range(adj_matrix.shape[0]):
        for jdx in range(adj_matrix.shape[0]):
            if adj_matrix[idx, jdx] == 1 and idx <= jdx:
                adj_list.append([idx, jdx])
    return adj_list

# Returns a list of successors by vertex, sensitive to the read strand
def to_adj_succ_list(adj_matrix, vertex_map, read_strand):
    succ_list = {}
    assert(adj_matrix.shape[0] == adj_matrix.shape[1])
    for idx in range(adj_matrix.shape[0]):
        succ_list[idx] = []
        for jdx in range(adj_matrix.shape[0]):
            if adj_matrix[idx, jdx] == 1:
                if read_strand == "+" and vertex_map[0][idx] <= vertex_map[0][jdx] or read_strand == "-" \
                   and vertex_map[1][idx] >= vertex_map[1][jdx]:
                    succ_list[idx].append(jdx)
    return succ_list

# Prints adjacency list representation of matrix
def print_adj_list(adj_list):
    print("EDGES: ")
    for edge in adj_list:
        print("{} -> {}".format(edge[0], edge[1]))



# Translate a DNA sequence encoding a peptide to amino-acid sequence via RNA
def translate_dna_to_peptide(dna_str):
    codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'
    }
    dna_str = dna_str.upper()
    aa_str = list("X"*(len(dna_str)/3))
    for idx in range(0, len(dna_str), 3):
        codon = dna_str[idx:idx+3]
        if len(codon) < 3:
            break
        if "N" in codon:
            aa_str[idx/3] = 'X'
        else:
            aa_str[idx/3] = codontable[codon]

    return "".join(aa_str)



# Returns true if there is a stop codon in the sequence. All codons that are fully in the
#  interval [start_coord<->end_coord] are checked.
# seq: Nucleotide sequence of vertex/CDS region
# start_coord: Read start coordinate
# stop_coord: Read stop coordinate
# strand: Read direction, one of {"+","-"}
def has_stop_codon_initial(seq, start_coord, end_coord, strand):
    if strand == "+":
        assert(start_coord <= end_coord)
        substr = seq[start_coord-1:end_coord]
    else:  # strand=="-"
        assert(start_coord >= end_coord)
        substr = complementary_seq(seq[end_coord-1:start_coord][::-1])
    for idx in range(0, len(substr)-2, 3):
        nuc_frame = substr[idx:idx+3]
        if nuc_frame.lower() in ["tag", "taa", "tga"]:
            return True
    return False



# Returns the to-from peptide frame depending on read strand and read_frame of the emitting vertex, returns
#  the stop codon in the second result.
# <emitting_frame>: Read frame of the emitting exon
# <strand>: Read strand
# <peptide_prop_seq_ref>: Reference sequence of the propagating exon
# <peptide_accept_seq_ref>: Reference sequence of the accepting exon
# <peptide_prop_seq_mut>: Mutated sequence of the propagating exon
# <peptide_accept_seq_mut>: Mutated sequence of the accepting exon
# <peptide_prop_coord>: Gene coordinate of the propagating exon
# <peptide_accept_coord>: Gene coordinate of the accepting exon
def cross_peptide_result(emitting_frame, strand, peptide_prop_seq_ref, peptide_accept_seq_ref,
                         peptide_prop_seq_mut, peptide_accept_seq_mut, peptide_prop_coord,
                         peptide_accept_coord):
    assert(len(peptide_prop_seq_mut) == len(peptide_prop_seq_ref))
    assert(len(peptide_accept_seq_mut) == len(peptide_accept_seq_ref))
    comp_read_frame = (3-emitting_frame) % 3
    prop_rest = (len(peptide_prop_seq_mut)-emitting_frame) % 3
    accept_rest = (len(peptide_accept_seq_mut)-comp_read_frame) % 3
    cor_accept_rest = -1*accept_rest if accept_rest > 0 else len(peptide_accept_seq_mut)
    if strand == "+":
        peptide_dna_str_mut = peptide_prop_seq_mut[prop_rest:]+peptide_accept_seq_mut[:cor_accept_rest]
        peptide_dna_str_ref = peptide_prop_seq_ref[prop_rest:]+peptide_accept_seq_ref[:cor_accept_rest]
        peptide_start_coord_v1 = peptide_prop_coord[0]+prop_rest
        peptide_stop_coord_v1 = peptide_prop_coord[1]
        peptide_start_coord_v2 = peptide_accept_coord[0]
        peptide_stop_coord_v2 = peptide_accept_coord[1]-accept_rest
    else:  # strand=="-"
        peptide_dna_str_mut = complementary_seq(peptide_prop_seq_mut[::-1][prop_rest:]+peptide_accept_seq_mut[::-1][:cor_accept_rest])
        peptide_dna_str_ref = complementary_seq(peptide_prop_seq_ref[::-1][prop_rest:]+peptide_accept_seq_ref[::-1][:cor_accept_rest])
        peptide_stop_coord_v1 = peptide_prop_coord[1]-prop_rest
        peptide_start_coord_v1 = peptide_prop_coord[0]
        peptide_stop_coord_v2 = peptide_accept_coord[1]
        peptide_start_coord_v2 = peptide_accept_coord[0]+accept_rest

    assert(len(peptide_dna_str_mut) == len(peptide_dna_str_ref))
    peptide_mut = translate_dna_to_peptide(peptide_dna_str_mut)
    peptide_ref = translate_dna_to_peptide(peptide_dna_str_ref)
    stop_codon_mut = False
    assert(len(peptide_dna_str_mut) % 3 == 0)

    # Detect pre/post junction stop codon in the sequence
    for idx in range(0, len(peptide_dna_str_mut), 3):
        nuc_frame = peptide_dna_str_mut[idx:idx+3]
        if nuc_frame.lower() in ["tag", "taa", "tga"]:
            stop_codon_mut = True
            break

    return (peptide_mut, peptide_ref, peptide_start_coord_v1+1, peptide_stop_coord_v1,
            peptide_start_coord_v2+1, peptide_stop_coord_v2, stop_codon_mut)

# Returns true if there is a stop codon found spanning the two sequences
# seq_prop: Propagating sequence
# seq_accept: Accepting sequence
# read_frame: Read frame of propagating vertex
# strand: Direction of read strand
def has_stop_codon_cross(seq_prop, seq_accept, read_frame, strand):
    if strand == "+":
        check_seq = seq_prop[-read_frame:]+seq_accept
    else:  # strand=="-"
        check_seq = seq_prop[::-1][-read_frame:]+seq_accept[::-1]
    for idx in range(0, len(check_seq)-2, 3):
        nuc_frame = check_seq[idx:idx+3]
        if nuc_frame.lower() in ["tag", "taa", "tga"]:
            return True
    return False


# Yields the complementary DNA sequence
# dna_seq: Input nucleotide sequence
def complementary_seq(dna_seq):
    comp_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    comp_dict_keys = comp_dict.keys()
    return "".join(map(lambda nuc: comp_dict[nuc] if nuc in comp_dict_keys else nuc, dna_seq))



# Encodes chromosome to same cn
def encode_chromosome(in_num):
    convert_dict = {23: "X", 24: "Y", 25: "MT"}
    return convert_dict[in_num] if in_num in convert_dict.keys() else str(in_num)
