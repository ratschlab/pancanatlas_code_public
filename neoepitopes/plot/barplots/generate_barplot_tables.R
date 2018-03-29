#
# generate data tables for peptides/sites bar plots
#
args = commandArgs(T)

base_dir = args[1]
input_file_asns = args[2]
input_file_snvs = args[3]
output_prefix = args[4]


#read in donors and cancer types
cancertype_file=paste0(base_dir, 'data/donor_to_cancertype.tcga.txt')
cancertypes = read.table(cancertype_file, stringsAsFactors = F, row.names=1)
colnames(cancertypes) = c('SET')

#main bar plot
data_barplot = paste0(output_prefix, '_expression-filtered.tsv')

ASNs = read.table(input_file_asns, row.names=1, stringsAsFactors=F)
colnames(ASNs) = c('ASN.unfiltered', 'ASN.filtered')

SNVs = read.table(input_file_snvs, row.names=1, stringsAsFactors=F)
colnames(SNVs) = c('SNV')

M = merge(ASNs, SNVs, by=0)
M = merge(M, cancertypes, by.x=1, by.y=0)
M$UNION.unfiltered = M$ASN.unfiltered + M$SNV
M$UNION.filtered = M$ASN.filtered + M$SNV
colnames(M)[1] = c('Sample')


D = data.frame(SET=character(0), TYPE=character(0), COUNT=numeric(0), stringsAsFactors=F)
ct = 'BRCA'
D[nrow(D)+1,] = c(ct, 'ASN', mean(M[which(M$SET == ct), 'ASN.filtered']))
D[nrow(D)+1,] = c(ct, 'SNV', mean(M[which(M$SET == ct), 'SNV']))
D[nrow(D)+1,] = c(ct, 'UNION', mean(M[which(M$SET == ct), 'UNION.filtered']))
ct = 'OV'
D[nrow(D)+1,] = c(ct, 'ASN', mean(M[which(M$SET == ct), 'ASN.filtered']))
D[nrow(D)+1,] = c(ct, 'SNV', mean(M[which(M$SET == ct), 'SNV']))
D[nrow(D)+1,] = c(ct, 'UNION', mean(M[which(M$SET == ct), 'UNION.filtered']))
ct = 'Total'
D[nrow(D)+1,] = c(ct, 'ASN', mean(M[,'ASN.filtered']))
D[nrow(D)+1,] = c(ct, 'SNV', mean(M[, 'SNV']))
D[nrow(D)+1,] = c(ct, 'UNION', mean(M[, 'UNION.filtered']))

write.table(D, file=data_barplot, quote=F, row.names=F, col.names=T, sep='\t')

#supplemental bar plot without expression filtering
data_barplot_S = paste0(output_prefix, '.tsv')

D = data.frame(SET=character(0), TYPE=character(0), COUNT=numeric(0), stringsAsFactors=F)
ct = 'BRCA'
D[nrow(D)+1,] = c(ct, 'ASN', mean(M[which(M$SET == ct), 'ASN.unfiltered']))
D[nrow(D)+1,] = c(ct, 'SNV', mean(M[which(M$SET == ct), 'SNV']))
D[nrow(D)+1,] = c(ct, 'UNION', mean(M[which(M$SET == ct), 'UNION.unfiltered']))
ct = 'OV'
D[nrow(D)+1,] = c(ct, 'ASN', mean(M[which(M$SET == ct), 'ASN.unfiltered']))
D[nrow(D)+1,] = c(ct, 'SNV', mean(M[which(M$SET == ct), 'SNV']))
D[nrow(D)+1,] = c(ct, 'UNION', mean(M[which(M$SET == ct), 'UNION.unfiltered']))
ct = 'Total'
D[nrow(D)+1,] = c(ct, 'ASN', mean(M[,'ASN.unfiltered']))
D[nrow(D)+1,] = c(ct, 'SNV', mean(M[, 'SNV']))
D[nrow(D)+1,] = c(ct, 'UNION', mean(M[, 'UNION.unfiltered']))

write.table(D, file=data_barplot_S, quote=F, row.names=F, col.names=T, sep='\t')

