''' In this script we output a table that collects the exons associated with each peptide'''

import glob
import sys
import os
import csv

base_dir = sys.argv[1]
tcga_base_path = sys.argv[2]
donor_file = sys.argv[3]
data_dir = sys.argv[4]

plot_dir=os.path.join(base_dir, "figures/")
exon_back_table_path=os.path.join(plot_dir,"data/data_recurrence_plot.tsv")

def main():
    all_donors=[] 
    with open(donor_file) as fid_in:
        for line in fid_in:
            all_donors.append(line.strip())

    print("Number of donors: {}".format(len(all_donors)))
    fp_out=open(exon_back_table_path,'w')
    csvfp_out=csv.writer(fp_out,delimiter="\t")
    csvfp_out.writerow(["DONOR", "PEPTIDE", "MHC_COUNTER", "SOURCE","TRANSCRIPT_SOURCE", "TYPE","MUTATE_MODE", "CHR", "GENE", "STRAND", 
                        "EXON1_START", "EXON1_END", "EXON2_START", "EXON2_END"])

    for donor in all_donors:
        jnct_file_ann=os.path.join(data_dir,"{}.gtex.unique_novel_junction_binders.2.0.annotated.new.txt".format(donor))
        n_jnct_peptides=0

        with open(jnct_file_ann,'r') as fp:
            for idx,line in enumerate(csv.reader(fp, delimiter="\t")):
                if idx==0:
                    continue
                n_jnct_peptides+=1
                peptide=line[0].strip()
                exon_spec=line[-1].strip()

                if exon_spec=="-1":
                    continue
                
                source_proteins=exon_spec.split(';')

                for prot in source_proteins:
                    if prot.split('_')[0][:4]=="TCGA" or prot.split('_')[0][:4]=="REFE":
                        exon_bounds=prot.split('_')[8:12]
                        chromosome=prot.split('_')[2].strip()
                        gene=prot.split("_")[1].strip()
                        strand=prot.split("_")[3].strip()
                        source=prot.split("_")[0].strip()
                        tf_source=prot.split("_")[5].strip()
                        mutate_mode=prot.split("_")[6].strip()
                        mhc_counter=-1
                    else:
                        exon_bounds=prot.split('_')[9:13]
                        chromosome=prot.split('_')[3].strip()
                        gene=prot.split("_")[2].strip()
                        strand=prot.split("_")[4].strip()
                        source=prot.split("_")[1].strip()
                        tf_source=prot.split("_")[6].strip()
                        mutate_mode=prot.split("_")[7].strip()
                        mhc_counter=int(prot.split("_")[0].strip())

                    assert(len(exon_bounds)==4)
                    csvfp_out.writerow([donor,peptide,mhc_counter,source,tf_source,"JNCT",mutate_mode,chromosome,gene,strand,
                                        int(exon_bounds[0]),int(exon_bounds[1]),int(exon_bounds[2]),int(exon_bounds[3])])

        print("Donor {}, JNCT: {}".format(donor,n_jnct_peptides))
    
    fp_out.close()

if __name__=="__main__":
    main()
