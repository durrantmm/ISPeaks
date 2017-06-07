Authors: Matt Durrant, Eli Moss, Ryan McKinney
Date: June 6 2017
Class: Statistical and Machine Learning Methods for Genomics

# Welcome to ISPeaks - The insertion sequence peak caller.

There are two folders of interest: `scripts and data`

`scripts` contains all of the most important code used for this project. The contents are briefly described as follows:

* `scripts`
  * `1.peak_shift.py` - This script calculates the appropriate distance to shift forward and reverse reads to maximize
the depth and improve our ability to identify peaks. Example usage:
   ```
  python scripts/1.peak_shift.py sample.in.bam cross_correlations.out.txt offset.txt corrected.out.bam
   ```
   
  * `2.significant_peaks.py` - This script calculates the signficant peaks within a given bam file. The output is in 
ENCODE narrowPeaks format. Example usage:
   ```
   python scripts/2.significant_peaks.py shifted.in.bam peaks.out.tsv peakshift_length 
   ```

   * `3.merge_peaks.py` - Takes two or more narrowPeaks files and merges the peaks, creating a union of
  all of the peaks. Example usage:
   ```
   python scripts/3.merge_peaks.py peaks.in.1.tsv peaks.in.2.tsv peaks.in.3.tsv ... > peaks.out.tsv 
   ```
 
   * `4.bootstrap_peaks.py` - This script performs bootstrap resampling of a bam file to produce fresh bootstrap
    realizations for downstream analysis. Useful for assessing signficance of differential peak abundants
    at a later step. Example usage:
   ```
   python scripts/4.bootstrap_peaks.py merged_peaks.in.tsv shifted.in.bam > bootstrap.read_counts.out.tsv 
   ```
  
   * `5.differential_peaks.R` - This script calculates the signficant peaks within a given bam file. The output is in 
ENCODE narrowPeaks format. Example usage:
   ```
   Rscript scripts/5.differential_peaks.R shifted.in.bam peaks.out.tsv peakshift_length 
   ```