library(cluster)
library(Biobase)
library(qvalue)
library(fastcluster)
options(stringsAsFactors = FALSE)
NO_REUSE = F

# try to reuse earlier-loaded data if possible
if (file.exists("1_GO-0005509~calcium_ion_binding.txt.RData") && ! NO_REUSE) {
    print('RESTORING DATA FROM EARLIER ANALYSIS')
    load("1_GO-0005509~calcium_ion_binding.txt.RData")
} else {
    print('Reading matrix file.')
    primary_data = read.table("1_GO-0005509~calcium_ion_binding.txt", header=T, com='', row.names=1, check.names=F, sep='\t')
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
write.table(data, file="1_GO-0005509~calcium_ion_binding.txt.log2.ZscaleRows.dat", quote=F, sep='	');
if (nrow(data) < 2) { stop("

**** Sorry, at least two rows are required for this matrix.

");}
if (ncol(data) < 2) { stop("

**** Sorry, at least two columns are required for this matrix.

");}
sample_dist = dist(t(data), method='euclidean')
hc_samples = NULL
gene_cor = NULL
gene_dist = dist(data, method='euclidean')
if (nrow(data) <= 1) { message('Too few genes to generate heatmap'); quit(status=0); }
hc_genes = hclust(gene_dist, method='complete')
heatmap_data = data
pdf("1_GO-0005509~calcium_ion_binding.txt.log2.ZscaleRows.genes_vs_samples_heatmap.pdf")
heatmap.3(heatmap_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=F, col=myheatcol, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=0.5, margins=c(10,10), cex.main=0.75, main=paste("samples vs. features
", "1_GO-0005509~calcium_ion_binding.txt.log2.ZscaleRows" ) )
dev.off()
