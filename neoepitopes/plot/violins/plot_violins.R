# plot kmer expression violins
library(ggplot2)

args = commandArgs(T)

base_dir = args[1]
plot_dir = paste0(base_dir, "/output/figures/")
sample_list = paste0(base_dir, "/donors.txt"
data_dir = paste0(base_dir,  "/output/cptac_mhc/cptac_validated_peptides/all/combined/"


# figure dimensions and colors
p_height = 6
p_width = 8

snv_color = "forestgreen"
asn_color = "blue"
bg_color = "darkgray"


modes = c('asn', 'snv', 'bg')
exclude_snv_set = TRUE #exclude SNV 9mers from ASN list

#
# read in data
#
samples = read.table(sample_list, stringsAsFactors=F)[,1]
D_full = NULL

for (mo in modes) {
    D = NULL
    Rdata_file = paste0(data_dir, "/donors.all.", mo, ".Rdata")
    print(mo)
    load(Rdata_file)

    # reduce to sample list
    D = subset(D, SAMPLE %in% samples)

    if (mo == "asn") {
        if (exclude_snv_set){
            D = subset(D, IN_SNV_SET == 0)
        } 
        D = subset(D, select = -c(IN_SNV_SET))
    }

    D$type = mo

    if (is.null(D_full)){
        D_full = D
    } else {
        D_full = rbind(D_full, D)
    }
}


D = D_full

#
# analyse data
#
print("Finished reading data ... Starting analysis.")

# exclude_0s: only consider 9-mers with expression > 0
for (exclude_0s in c(TRUE)){
    print(exclude_0s)
    zeros = ""
    if (exclude_0s){
        D = subset(D_full, SCORE > 0)
        zeros = ".exclude_0s"
    } else {
        D = D_full
        zeros = ""
    } 
       
    all_lists = list(Dmo = list(asn=NULL, snv=NULL), D_noGTEX = list(asn=NULL, snv=NULL))

    for (mo in modes){
        Dmo = subset(D, type == mo)
        print(paste0(mo, " Dmo ",dim(Dmo), " cptac ", sum(Dmo$IN_CPTAC)))
        
        D_noGTEX = subset(Dmo, IN_GTEX == 0)
        print(paste0("D_noGTEX ",dim(D_noGTEX), " cptac ", sum(D_noGTEX$IN_CPTAC)))

        all_lists[["Dmo"]][[mo]] = Dmo
        all_lists[["D_noGTEX"]][[mo]] = D_noGTEX
    } 

    for (L in c("D_noGTEX")){         
        print(L)
        D_list = all_lists[[L]]

        K = D_list[["asn"]][,c("SCORE", "IN_CPTAC")]
        x = dim(K)[1]
        K = rbind(K, D_list[["snv"]][, c("SCORE",  "IN_CPTAC")])
        y = dim(K)[1]
        K = rbind(K, D_list[['bg']][, c("SCORE",  "IN_CPTAC")])
        K$SET = "OVERALL"
        K$SET[1:y] = "SNV"
        K$SET[1:x] = "ASN"
        K$SET = as.factor(K$SET)
        K$SET = factor(K$SET,levels=c("ASN", "SNV", "OVERALL"))
        print(table(K$SET))
        K$LOG_SCORE = log10(K$SCORE + 0.0000001)

        pdf(paste0(plot_dir, '/', L, '.score_distribution.log_exp.2violins', zeros, '.test.pdf'), height=p_height, width=p_height)
        p = ggplot(K, aes(x=SET, y=LOG_SCORE, color=SET)) + 
        geom_violin(trim=FALSE, fill='gray', linetype='dotted') + 
        geom_violin(data=K[which(K$IN_CPTAC == 1),], mapping=aes(x=SET, y=LOG_SCORE, color=SET), trim=FALSE, fill=NA) + 
        scale_color_manual(values = c(asn_color, snv_color, bg_color)) +
        xlab("Set") + 
        ylab("log10(expression)") +
        theme_bw() + 
        theme(axis.text=element_text(size=rel(2)), axis.title=element_text(size=rel(2))) + 
        theme(legend.position="none") 
        plot(p)
        dev.off()
    }
}
