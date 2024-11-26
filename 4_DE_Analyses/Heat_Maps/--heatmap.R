library(cluster)
library(Biobase)
library(qvalue)
library(fastcluster)
options(stringsAsFactors = FALSE)
NO_REUSE = F

# try to reuse earlier-loaded data if possible
if (file.exists("--heatmap.RData") && ! NO_REUSE) {
    print('RESTORING DATA FROM EARLIER ANALYSIS')
    load("--heatmap.RData")
} else {
    print('Reading matrix file.')
    primary_data = read.table("--heatmap", header=T, com='', row.names=1, check.names=F, sep='\t')
    primary_data = as.matrix(primary_data)
}
source("/Applications/trinityrnaseq-v2.15.2/Analysis/DifferentialExpression/R/heatmap.3.R")
source("/Applications/trinityrnaseq-v2.15.2/Analysis/DifferentialExpression/R/misc_rnaseq_funcs.R")
source("/Applications/trinityrnaseq-v2.15.2/Analysis/DifferentialExpression/R/pairs3.R")
source("/Applications/trinityrnaseq-v2.15.2/Analysis/DifferentialExpression/R/vioplot2.R")
data = primary_data
myheatcol = colorpanel(75, 'blue','black','yellow')
data = data[, c('C_38_FPKM','C_51_FPKM','E_38_FPKM','E_51_FPKM'), drop=F ]
sample_types = colnames(data)
nsamples = length(sample_types)
sample_colors = rainbow(nsamples)
sample_type_list = list()
for (i in 1:nsamples) {
    sample_type_list[[sample_types[i]]] = sample_types[i]
}
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
initial_matrix = data # store before doing various data transformations
data = log2(data+1)
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
sampleAnnotations = matrix(ncol=ncol(data),nrow=nsamples)
for (i in 1:nsamples) {
  sampleAnnotations[i,] = colnames(data) %in% sample_type_list[[sample_types[i]]]
}
sampleAnnotations = apply(sampleAnnotations, 1:2, function(x) as.logical(x))
sampleAnnotations = sample_matrix_to_color_assignments(sampleAnnotations, col=sample_colors)
rownames(sampleAnnotations) = as.vector(sample_types)
colnames(sampleAnnotations) = colnames(data)
data = as.matrix(data) # convert to matrix
# Z-scale the genes across all the samples for PCA
zscaled_data = t(scale(t(data), scale=T))
data = zscaled_data
data = na.omit(data)
write.table(data, file="--heatmap.log2.ZscaleRows.dat", quote=F, sep='	');
gene_cor = NULL
