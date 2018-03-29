''' Generates hdf5 file from pre_data tsv for recurrence plot'''

import os
import os.path
import glob
import csv
import numpy as np
import pandas as pd


base_dir = sys.argv[1]

plot_dir=os.path.join(base_dir, "figures/")
pre_data=os.path.join(plot_dir,"data/data_recurrence_plot.tsv")
donor_to_cancertype_file = os.path.join(base_dir, 'data/donor_to_cancertype.tcga.txt')
output_file = os.path.join(base_dir, "output/figures/data/data_recurrence_plot.h5")


def parse_donor_cancer_type_map():
    donor_cancer_dict={}
    with open(donor_to_cancertype_file,'r') as fp:
        csvfp=csv.reader(fp,delimiter='\t')

        for row in csvfp:
            donor_id=row[0].strip()
            cancer_type=row[1].strip()
            assert(donor_id not in donor_cancer_dict)
            donor_cancer_dict[donor_id]=cancer_type

    return donor_cancer_dict

def main():
    all_kmers=[]
    all_donor_ids=[]
    donor_kmer_map={}
    kmer_donor_map={}
    donor_cancer_map=parse_donor_cancer_type_map()
    cancer_list=list(sorted(set(donor_cancer_map.values())))

    with open(pre_data,'r') as fp:
        tsvfp=csv.reader(fp,delimiter='\t')
        idx=0
        next(tsvfp, None)  # skip the headers
        for row in tsvfp:
            idx+=1
            donor_id=row[0].strip()
            kmer=row[1].strip()
            all_kmers.append(kmer)
            all_donor_ids.append(donor_id)

            if donor_id not in donor_kmer_map:
                donor_kmer_map[donor_id]=[kmer]
            else:
                donor_kmer_map[donor_id].append(kmer)

            if kmer not in kmer_donor_map:
                kmer_donor_map[kmer]=[donor_id]
            else:
                kmer_donor_map[kmer].append(donor_id)

    # Convert everything to sets, we are only interested if the k-mer occurs at least once
    for kmer in kmer_donor_map:
        kmer_donor_map[kmer]=set(kmer_donor_map[kmer])
    for donor in donor_kmer_map:
        donor_kmer_map[donor]=set(donor_kmer_map[donor])

    all_kmers=sorted(list(set(all_kmers)))
    print("Number of all kmers: {}".format(len(all_kmers)))
    recurrent_kmers=[]
    
    for kmer in all_kmers:
        if len(kmer_donor_map[kmer])>1:
            recurrent_kmers.append(kmer)

    recurrent_kmers=sorted(recurrent_kmers,key=lambda elem: len(kmer_donor_map[elem]))
    print("Number of recurrent K-mers: {}".format(len(recurrent_kmers)))
    all_donor_ids=sorted(list(set(all_donor_ids)))
    recurrent_donors=[]

    for donor_id in all_donor_ids:
        donor_kmers=set(donor_kmer_map[donor_id])
        for cont_kmer in donor_kmers:
            if cont_kmer in recurrent_kmers:
                recurrent_donors.append(donor_id)
                break

    # Update donor-kmer map such that only recurrent k-mers are taken into account
    for donor in donor_kmer_map:
        donor_kmer_map[donor] = list(filter(lambda elem: elem in recurrent_kmers, donor_kmer_map[donor]))

    recurrent_donors=sorted(recurrent_donors, key=lambda elem: (cancer_list.index(donor_cancer_map[elem]), len(donor_kmer_map[elem])))
    print("Number of unique donors: {}".format(len(all_donor_ids)))
    print("Number of recurrent donors: {}".format(len(recurrent_donors)))

    data_arr_dict={}

    kmer_idx_map={}
    for kmer_idx,kmer in enumerate(recurrent_kmers):
        kmer_idx_map[kmer]=kmer_idx

    donor_idx_map={}
    for donor_idx,donor in enumerate(recurrent_donors):
        donor_idx_map[donor]=donor_idx

    data_arr_matrix=np.zeros((len(recurrent_kmers), len(recurrent_donors)), dtype=np.int)
    data_norm_matrix=np.zeros((len(recurrent_kmers), len(recurrent_donors)), dtype=np.int)

    for donor_id in recurrent_donors:
        for kmer in donor_kmer_map[donor_id]:
            data_arr_matrix[kmer_idx_map[kmer], donor_idx_map[donor_id]]=1+cancer_list.index(donor_cancer_map[donor_id])
            data_norm_matrix[kmer_idx_map[kmer], donor_idx_map[donor_id]]=1

    print("Collected donor data...")
    df_nonplot=pd.DataFrame(data=data_norm_matrix,columns=recurrent_donors,index=recurrent_kmers)
    df_nonplot.to_hdf(output_file,'data')


if __name__=="__main__":
    main()
