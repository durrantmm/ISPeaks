#!/usr/bin/env Rscript

require(boot)
#require(DESeq2)
#require(ggplot2)
#require("pheatmap")

#args <- commandArgs(trailingOnly = TRUE)
#file1 = args[1]
#file2 = args[2]

#for dev
setwd('~/scg4_moss/projects/coursework/gene_245/ISPeaks/')
file1 = 'data/peaks/6753_12-15-15.47678.IS614.peaks.readcounts.tsv'
file2 = 'data/peaks/6753_3-1-16.47678.IS614.peaks.readcounts.tsv'
d.1 = read.table(file1, headers = T)
d.2 = read.table(file2, headers = T)

#bootstrap a number of samples from each

#setup to call Deseq2
group = c(1, 2)

#some day 20 samples are controls!
day[8:9] = '20ctrl'
#include them in with the day 0 timepoints
day[day == '20ctrl'] = 0
day[day == '14'] = 13

id = colnames(d)
metadata = data.frame(cbind(day, id))


input = d
data = DESeqDataSetFromMatrix(countData = input, colData = metadata, design = ~ day)
proc = DESeq(data, test = "Wald")
res = results(proc)
