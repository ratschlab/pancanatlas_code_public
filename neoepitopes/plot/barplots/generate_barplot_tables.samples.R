#
# generate data tables for sample fraction bar plots
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
colnames(ASNs) = c('unfiltered', 'filtered')
ASNs$unfiltered[which(ASNs$unfiltered > 0)] = 1
ASNs$filtered[which(ASNs$filtered > 0)] = 1

SNVs = read.table(input_file_snvs, row.names=1, stringsAsFactors=F)
colnames(SNVs) = c('SNV')
SNVs = merge(SNVs, cancertypes, by=0)
SNVs$SNV[which(SNVs$SNV > 0)] = 1

ASNs = read.table(input_file_asns, row.names=1, stringsAsFactors=F)
colnames(ASNs) = c('ASN.unfiltered', 'ASN.filtered')

SNVs = read.table(input_file_snvs, row.names=1, stringsAsFactors=F)
colnames(SNVs) = c('SNV')

M = merge(ASNs, SNVs, by=0)
M = merge(M, cancertypes, by.x=1, by.y=0)
M$ASN.unfiltered[which(M$ASN.unfiltered > 0)] = 1
M$ASN.filtered[which(M$ASN.filtered > 0)] = 1
M$SNV[which(M$SNV > 0)] = 1
M$UNION.unfiltered = 1 * (M$ASN.unfiltered | M$SNV)
M$UNION.filtered = 1 * (M$ASN.filtered | M$SNV)
colnames(M)[1] = c('Sample')



D = data.frame(SET=character(0), TYPE=character(0), COUNT=numeric(0), stringsAsFactors=F)
ct = 'BRCA'
nof_samples = length(which(M$SET == ct))
D[nrow(D)+1,] = c(ct, 'ASN', sum(M[which(M$SET == ct), 'ASN.filtered'])/nof_samples)
D[nrow(D)+1,] = c(ct, 'SNV', sum(M[which(M$SET == ct), 'SNV'])/nof_samples)
D[nrow(D)+1,] = c(ct, 'UNION', sum(M[which(M$SET == ct), 'UNION.filtered'])/nof_samples)
ct = 'OV'
nof_samples = length(which(M$SET == ct))
D[nrow(D)+1,] = c(ct, 'ASN', sum(M[which(M$SET == ct), 'ASN.filtered'])/nof_samples)
D[nrow(D)+1,] = c(ct, 'SNV', sum(M[which(M$SET == ct), 'SNV'])/nof_samples)
D[nrow(D)+1,] = c(ct, 'UNION', sum(M[which(M$SET == ct), 'UNION.filtered'])/nof_samples)
ct = 'Total'
nof_samples = length(M$SET)
D[nrow(D)+1,] = c(ct, 'ASN', sum(M[,'ASN.filtered'])/nof_samples)
D[nrow(D)+1,] = c(ct, 'SNV', sum(M[, 'SNV'])/nof_samples)
D[nrow(D)+1,] = c(ct, 'UNION', sum(M[, 'UNION.filtered'])/nof_samples)

write.table(D, file=data_barplot, quote=F, row.names=F, col.names=T, sep='\t')

#supplemental bar plot without expression filtering
data_barplot_S = paste0(output_prefix, '.tsv')

D = data.frame(SET=character(0), TYPE=character(0), COUNT=numeric(0), stringsAsFactors=F)
ct = 'BRCA'
nof_samples = length(which(M$SET == ct))
D[nrow(D)+1,] = c(ct, 'ASN', sum(M[which(M$SET == ct), 'ASN.unfiltered'])/nof_samples)
D[nrow(D)+1,] = c(ct, 'SNV', sum(M[which(M$SET == ct), 'SNV'])/nof_samples)
D[nrow(D)+1,] = c(ct, 'UNION', sum(M[which(M$SET == ct), 'UNION.unfiltered'])/nof_samples)
ct = 'OV'
nof_samples = length(which(M$SET == ct))
D[nrow(D)+1,] = c(ct, 'ASN', sum(M[which(M$SET == ct), 'ASN.unfiltered'])/nof_samples)
D[nrow(D)+1,] = c(ct, 'SNV', sum(M[which(M$SET == ct), 'SNV'])/nof_samples)
D[nrow(D)+1,] = c(ct, 'UNION', sum(M[which(M$SET == ct), 'UNION.unfiltered'])/nof_samples)
ct = 'Total'
nof_samples = length(M$SET)
D[nrow(D)+1,] = c(ct, 'ASN', sum(M[,'ASN.unfiltered'])/nof_samples)
D[nrow(D)+1,] = c(ct, 'SNV', sum(M[, 'SNV'])/nof_samples)
D[nrow(D)+1,] = c(ct, 'UNION', sum(M[, 'UNION.unfiltered'])/nof_samples)

write.table(D, file=data_barplot_S, quote=F, row.names=F, col.names=T, sep='\t')
