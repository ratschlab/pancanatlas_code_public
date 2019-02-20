
''' IO-related utility functions'''

from __future__ import print_function

import csv
import datetime
import os

# Saves a list to CSV file, one element per line
def save_list_csv(in_lst, out_path):
    with open(out_path, 'w') as fp:
        for elem in in_lst:
            print(elem, file=fp)


# Setup local paths for writing the output FASTA file
def generate_output_fp(donor_list_name, donor_id, output_type, configs):
    output_base_dir=configs["output_path"]
    timestamp = datetime.datetime.now().strftime("%y%m%d")
    output_path = os.path.join(output_base_dir, "peptides_{}_{}".format(donor_list_name, timestamp))

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    output_fname = "{}_{}.fa".format(donor_id, output_type)
    output_path = os.path.join(output_path, output_fname)
    output_fp = open(output_path, "w")
    return output_fp

# Returns donor list to be processed
def donor_list(configs):
    donor_filepath=configs["donor_list_path"]
    donor_ids = []
    with open(donor_filepath, "r") as fp:
        csv_fp = csv.reader(fp, delimiter="\t")
        for row in csv_fp:
            first_elem = row[0].strip()

            # Chop off prefix if not in standard-format
            if not first_elem[:4] == "TCGA":
                donor_id = "-".join(first_elem.split("-")[1:])
            else:
                donor_id = first_elem

            donor_ids.append(donor_id)
    return sorted(donor_ids)

