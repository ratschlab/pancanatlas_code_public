''' Misc. utility functions'''

# Returns header labels corresponding to donor_id and mutation_mode
def header_labels(donor_id, mutation_mode):
    if mutation_mode is None:
        mutation_type = "REFERENCE"
    elif mutation_mode == "both":
        mutation_type = "GERM_SOMATIC"
    elif mutation_mode == "somatic_only":
        mutation_type = "SOMATIC"
    elif mutation_mode == "germline_only":
        mutation_type = "GERM"

    peptide_type = "REFERENCE" if donor_id == "ref" else donor_id
    return (peptide_type, mutation_type)
