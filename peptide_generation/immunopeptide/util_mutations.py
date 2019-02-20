''' Utility functions related to parsing of mutations file and applying them'''

import csv

import h5py
import numpy as np

from immunopeptide.util_genomic import encode_chromosome


# Applies mutations to the entire genome, given a list of mutations to consider
# seq_dict:
# apply_somatic: [TRUE,FALSE] should the somatic mutatione be applied?
# som_info: Look-up table for somatic mutations
# apply_germline: [TRUE,FALSE] should the germline mutations be applied?
# germ_info: Look-up table for germline mutations
def apply_mutations(seq_dict, apply_somatic=True, apply_germline=True, som_info=None, germ_info=None):
    personalized_seq_dict = {}

    # Convert the sequences to list of chars for random access
    for k, v in seq_dict.iteritems():
        personalized_seq_dict[k] = list(v)

    # Apply germline mutations to input if requested
    if apply_germline:

        # Apply each germline mutation
        for gl_mut in germ_info:
            mut_start = int(gl_mut["start_coord"])-1
            mut_stop = int(gl_mut["stop_coord"])-1
            mut_chr = gl_mut["chromosome"]
            assert(mut_start == mut_stop)
            ref_seq = seq_dict[mut_chr]
            assert(mut_start < len(ref_seq))
            mut_ref_nuc = gl_mut["ref_nuc"]
            assert(len(mut_ref_nuc) == 1)
            mut_new_nuc = gl_mut["mut_nuc"]
            assert(len(mut_new_nuc) == 1)
            mut_seq = personalized_seq_dict[mut_chr]
            assert(mut_ref_nuc == ref_seq[mut_start])

            # If an uncertainty placeholder exists, use the reference..
            if "*" not in mut_new_nuc:
                mut_seq[mut_start] = mut_new_nuc
                assert(personalized_seq_dict[mut_chr][mut_start] == mut_new_nuc)
                assert(seq_dict[mut_chr][mut_start] == mut_ref_nuc)

    # Apply somatic mutations to input if requested, after germline mutations
    if apply_somatic:

        # Apply each somatic mutation
        for som_mut in som_info:
            mut_start = int(som_mut["start_coord"])-1
            mut_stop = int(som_mut["stop_coord"])-1
            mut_chr = som_mut["chromosome"]
            mut_variant = som_mut["variant_type"]
            assert(mut_variant == "SNP")
            assert(mut_start == mut_stop)
            ref_seq = seq_dict[mut_chr]
            assert(mut_start < len(ref_seq))
            mut_ref_nuc = som_mut["ref_nuc"]
            assert(len(mut_ref_nuc) == 1)
            mut_new_nuc1 = som_mut["mut_nuc1"]
            mut_new_nuc2 = som_mut["mut_nuc2"]
            # Apply on the mutated allele, and if both or none are mutated take the
            # the first allele arbitrarily
            if mut_new_nuc1 != mut_ref_nuc:
                mut_new_nuc = mut_new_nuc1
            elif mut_new_nuc2 != mut_ref_nuc:
                mut_new_nuc = mut_new_nuc2
            else:
                mut_new_nuc = mut_new_nuc1
            assert(len(mut_new_nuc) == 1)
            mut_seq = personalized_seq_dict[mut_chr]
            assert(mut_ref_nuc == ref_seq[mut_start])
            mut_seq[mut_start] = mut_new_nuc
            assert(personalized_seq_dict[mut_chr][mut_start] == mut_new_nuc)
            assert(seq_dict[mut_chr][mut_start] == mut_ref_nuc)

    # Write mutated sequences to the output
    for k in personalized_seq_dict.keys():
        personalized_seq_dict[k] = "".join(personalized_seq_dict[k])

    return personalized_seq_dict




# Parses the somatic mutation file and builds a lookup-table with only the required information
def parse_somatic_mutation_info(donor_list, configs):
    som_path=configs["somatic_mutation_path"]
    somatic_table = {}
    somatic_pos_to_mutation_desc = {}

    with open(som_path,'r') as fp:
        csv_fp = csv.reader(fp, delimiter="\t")
        # Check if current row describes a donor
        for idx, row in enumerate(csv_fp):
            donor_id = "-".join(row[15].split("-")[0:3]).strip()
            if donor_id not in somatic_table:
                somatic_table[donor_id] = []
            variant_type = row[9].strip()
            #print("Donor-ID: {}".format(donor_id))
            #print("Variant type: {}".format(variant_type))
            if donor_id in donor_list and variant_type == "SNP":
                somatic_dict = {}
                somatic_dict["chromosome"] = row[4].strip()
                somatic_dict["start_coord"] = row[5].strip()
                somatic_dict["stop_coord"] = row[6].strip()
                somatic_dict["strand"] = row[7].strip()
                somatic_dict["variant_type"] = row[9].strip()
                somatic_dict["ref_nuc"] = row[10].strip()
                somatic_dict["mut_nuc1"] = row[11].strip()
                somatic_dict["mut_nuc2"] = row[12].strip()
                mut_desc = {}
                mut_desc["t_depth"] = row[39].strip()
                mut_desc["t_ref_count"] = row[40].strip()
                mut_desc["t_alt_count"] = row[41].strip()
                mut_desc["n_depth"] = row[42].strip()
                mut_desc["n_ref_count"] = row[43].strip()
                mut_desc["n_alt_count"] = row[44].strip()
                somatic_table[donor_id].append(somatic_dict)
                somatic_pos_to_mutation_desc[donor_id, int(somatic_dict["start_coord"])] = mut_desc
    return (somatic_table, somatic_pos_to_mutation_desc)



# Parses the germline mutation file and builds a look-up table with only the required information
def parse_germline_mutation_info(donor_list, configs):
    germ_path=configs["germline_mutation_path"]
    germline_table = {}

    # Load mutations from HDF5 file
    with h5py.File(germ_path, 'r') as h5f:
        allele_alt = np.array(h5f["allele_alt"])
        allele_ref = np.array(h5f["allele_ref"])
        allele_pos = np.array(h5f["pos"])
        gt_id = np.array(h5f["gtid"])
        sel_lst = []
        run_idx_map = {}
        run_idx = 0
        for idx, elem in np.ndenumerate(gt_id):
            tcga_key = '-'.join(gt_id[idx[0]].split('-')[0:3])

            if tcga_key in donor_list:
                sel_lst.append(idx[0])
                run_idx_map[tcga_key] = run_idx
                run_idx += 1

        gt_arr = np.array(h5f["gt"][:, sel_lst])

    # save_list_csv(sorted(run_idx_map.keys()),os.path.join(germ_path,"donors_remaining_170821_with_germline_info.txt"))
    assert(sorted(run_idx_map.keys()) == sorted(donor_list))

    # Build look-up dictionary...
    for donor in donor_list:
        donor_col_idx = run_idx_map[donor]

        for row_idx in np.arange(gt_arr.shape[0]):
            mut_status = gt_arr[row_idx, donor_col_idx]
            if mut_status <= 1.0:
                germline_dict = {}
                germline_dict["chromosome"] = encode_chromosome(int(allele_pos[row_idx, 0]))
                germline_dict["start_coord"] = allele_pos[row_idx, 1]
                germline_dict["stop_coord"] = allele_pos[row_idx, 1]
                germline_dict["ref_nuc"] = allele_ref[row_idx]
                germline_dict["mut_nuc"] = allele_alt[row_idx]

                if donor not in germline_table:
                    germline_table[donor] = []

                germline_table[donor].append(germline_dict)
    return germline_table




