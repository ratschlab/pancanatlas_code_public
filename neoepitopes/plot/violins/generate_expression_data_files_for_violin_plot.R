# generate expression input data for violin plot
# one file per mode: ASN, SNV, OVERALL (bg)
library(plyr)

args = commandArgs(T)

base_dir = args[1]
plot_dir = paste0(base_dir, "/output/figures/")
sample_list = paste0(base_dir, "/donors.txt"
data_dir = paste0(base_dir,  "/output/cptac_mhc/cptac_validated_peptides/all/combined/"

samples = read.table(sample_list, stringsAsFactors=F)[,1]


modes = c('asn', 'snv', 'bg')
for (mo in modes) {
    D = NULL
    Rdata_file = paste0(data_dir, "/donors.all.", mo, ".Rdata")
    print(mo)

    if (file.exists(Rdata_file)){
        print(paste0("Rdata file ", Rdata_file, " exists. Skipping ... "))
        next
    } else {
        for (sample in samples){
            Dchunk = NULL
            Dtmp = NULL
            sample_input_dir = paste0(input_dir, "/", mo, "s/")
            if (mo == "asn"){
                sample_Rdata_file = paste0(input_dir, "/asns/", sample, ".asn.all.Rdata")
                print(paste0("Reading ", sample_Rdata_file)) 
                tmp_env = new.env()
                load(sample_Rdata_file, tmp_env)
                X = get(ls(tmp_env), envir=tmp_env)
                rm(tmp_env)
            } else { # snv or bg
                sample_input_dir = paste0(input_dir, "/", mo, "s/")
                input = paste0(sample_input_dir, "/", sample, ".snv.table.txt")
                print(paste0("Reading ", input))
                Dchunk = NULL
                if (grepl("gz$", input)){
                    Dchunk = read.table(gzfile(input), header=T, stringsAsFactors=F, sep='\t')
                } else {
                    Dchunk = read.table(input, header=T, stringsAsFactors=F, sep='\t')
                }

                if (dim(Dchunk)[1] == 0){
                    print(paste0("Skipping empty chunk "))
                    next
                }

                Dchunk$SAMPLE=sample
                if (mo == "snv"){
                    #calculate norm expression
                    Dchunk$NORMALIZED_SNV_EXPRESSION = Dchunk$SNV_EXPRESSION/Dchunk$LIBRARY_SIZE
                    x = which(colnames(Dchunk) == "NORMALIZED_SNV_EXPRESSION")
                    colnames(Dchunk)[x] = "SCORE"
                } else { #bg
                    #calculate norm expression
                    Dchunk$NORMALIZED_EXPRESSION = Dchunk$EXPRESSION/Dchunk$LIBRARY_SIZE
                    x = which(colnames(Dchunk) == "NORMALIZED_EXPRESSION")
                    colnames(Dchunk)[x] = "SCORE"
                    Dchunk$IN_CPTAC = 0 # no CPTAC curve for bg/overall in violin plot
                }

                E = Dchunk[,c('KMER_SEQ', 'BIND_STRENGTH', 'SCORE', 'IN_GTEX', 'IN_CPTAC', 'SAMPLE')]
                X = ddply(E, .(KMER_SEQ), function(x)x[which.max(x$SCORE)[1],])

                Dtmp = X
            }
       

            if (is.null(dim(D))){
                D = X
            } else {
                D = rbind(D, X)
            }
        }

        # save data frame
        save(D, file=Rdata_file)
    }
}
