setwd('~/scg4_moss/projects/coursework/gene_245/ISPeaks/')

for (file in list.files('output/peakshift_outputs', pattern = '*.txt')){
	print(file)
	pdf(paste('images/peakshift_correlations/', file, '.pdf', sep = ''), height = 4, width = 3)
	cors = scan(paste('output/peakshift_outputs/', file, sep = ''))
	cors = cors[(length(cors)/2 - 1000):(length(cors)/2 )]
	plot(cors, type = 'l', xlab = 'Offset (Forward and Reverse)', ylab = 'Forward and Reverse Depth Cross-Correlation')
	abline(v = which.max(cors), lty = 3, col = 'blue')
	print(which.max(cors))
	dev.off()
}
