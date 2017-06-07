#!/usr/bin/env Rscript

#install deseq2
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

require(DESeq2)
#require(ggplot2)
require("pheatmap")
require(reshape2)

#args <- commandArgs(trailingOnly = TRUE)
#file1 = args[1]
#file2 = args[2]

#for dev
setwd('~/scg4_moss/projects/coursework/gene_245/ISPeaks/')
file1 = 'FinalCode/output/4.bootstrapped_counts/sim1.counts.tsv'
file2 = 'FinalCode/output/4.bootstrapped_counts/sim2.counts.tsv'
#file = 'output/bootstrap_resamples/6753_bootstrapped_counts.tsv'
#file = 'output/bootstrap_resamples/diff_sim.merged.bootstraps.tsv'
d.1 = read.table(file1, header = T)
d.1 = cbind(rep('samp1', nrow(d.1)), d.1)
colnames(d.1)[1] = "Sample"
d.2 = read.table(file2, header = T)
d.2 = cbind(rep('samp2', nrow(d.1)), d.2)
colnames(d.2)[1] = "Sample"
d.long = rbind(d.1, d.2)
#d.long = d.long[d.long$Sampling %in% c('EMPIRICAL', 'BOOTSTRAP1', 'BOOTSTRAP2', 'BOOTSTRAP3'),]
d.long.2vars = data.frame(paste(d.long$Sample, d.long$Sampling, sep = '_'), paste(d.long$Contig, d.long$PeakStart, sep = '_'), d.long$ReadCount)
colnames(d.long.2vars) = c('Sample', 'Peak', 'ReadCount')
d.wide = reshape(d.long.2vars, idvar = "Sample", timevar = "Peak", direction = "wide", new.row.names = unique(paste(d.long$Sample, d.long$Sampling, sep = '_')))
d.wide = d.wide[,2:ncol(d.wide)]
#setup to call Deseq2
group = sapply(rownames(d.wide), function(x) strsplit(x, '_')[[1]][1])
id = rownames(d.wide)
metadata = data.frame(cbind(group, id))

input = t(d.wide)
data = DESeqDataSetFromMatrix(countData = input, colData = metadata, design = ~ group)
proc = DESeq(data, test = "Wald")
res = results(proc)

betas = coef(proc)
topGenes <- head(order(res$padj),100)
mat <- betas[topGenes, -1]#[,c(1,4,5,6,2,3)]
thr <- 3
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
rownames(mat) = sapply(rownames(mat), function(x) gsub('__', '  ', x))
rownames(mat) = sapply(rownames(mat), function(x) strsplit(x, ' ')[[1]][1])
#pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),cluster_col=FALSE)

#jpeg('images/simulated deseq2 longitudinal.jpeg', width = 4, height = 4, unit = 'in', res = 160)
#plot(-log10(res$pvalue/(0.05/nrow(d.wide))), ylab = "-log10 Adjusted p-value", xlab = 'Peaks Sorted by Position', main = "Simulated Longitudinal Data")
plot(-log10(res$pvalue/(0.05/nrow(d.wide))), ylab = "-log10 Adjusted p-value", xlab = 'Peaks Sorted by Position', main = "Clinical Longitudinal Data")
#dev.off()

heatmap(as.matrix(d.wide))

plot(-log10(res$pvalue))
abline(h = -log10(0.05/nrow(d.wide)), lty = 3, col = 'darkgreen')
