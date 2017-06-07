#!/usr/bin/env bash

# Executing all of these cod

# Making the output folders
mkdir output
mkdir output/1.bams_shifted
mkdir output/2.significant_peaks
mkdir output/3.merged_peaks
mkdir output/4.bootstrapped_counts

# Shift the bam files
python scripts/1.peak_shift.py data/sim1.bacteroides_caccae.bam output/1.bams_shifted/sim1.bacteroides_caccae.cross_correl.txt output/1.bams_shifted/sim1.bacteroides_caccae.offset.txt output/1.bams_shifted/sim1.bacteroides_caccae.shifted.bam
python scripts/1.peak_shift.py data/sim1.bacteroides_caccae.bam output/1.bams_shifted/sim2.bacteroides_caccae.cross_correl.txt output/1.bams_shifted/sim2.bacteroides_caccae.offset.txt output/1.bams_shifted/sim2.bacteroides_caccae.shifted.bam

# Sort the sam files
samtools sort output/1.bams_shifted/sim1.bacteroides_caccae.shifted.bam > output/1.bams_shifted/sim1.bacteroides_caccae.shifted.sorted.bam
samtools sort output/1.bams_shifted/sim2.bacteroides_caccae.shifted.bam > output/1.bams_shifted/sim2.bacteroides_caccae.shifted.sorted.bam

# Index the sam files
samtools index output/1.bams_shifted/sim1.bacteroides_caccae.shifted.sorted.bam
samtools index output/1.bams_shifted/sim2.bacteroides_caccae.shifted.sorted.bam

# Call the peaks
python scripts/2.significant_peaks.py output/1.bams_shifted/sim1.bacteroides_caccae.shifted.sorted.bam output/2.significant_peaks/sim1.bacteroides_caccae.peaks.tsv 250
python scripts/2.significant_peaks.py output/1.bams_shifted/sim2.bacteroides_caccae.shifted.sorted.bam output/2.significant_peaks/sim2.bacteroides_caccae.peaks.tsv 250

# Merge the peaks
python scripts/3.merge_peaks.py output/2.significant_peaks/sim1.bacteroides_caccae.peaks.tsv output/2.significant_peaks/sim2.bacteroides_caccae.peaks.tsv >  output/3.merged_peaks/sim.merged.peaks.tsv

# Bootstrap the data to get bootstrapped read counts
python scripts/4.bootstrap_peaks.py output/3.merged_peaks/sim.merged.peaks.tsv output/1.bams_shifted/sim1.bacteroides_caccae.shifted.sorted.bam > output/4.bootstrapped_counts/sim1.counts.tsv
python scripts/4.bootstrap_peaks.py output/3.merged_peaks/sim.merged.peaks.tsv output/1.bams_shifted/sim2.bacteroides_caccae.shifted.sorted.bam > output/4.bootstrapped_counts/sim2.counts.tsv