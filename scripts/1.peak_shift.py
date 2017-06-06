#!/usr/bin/env python

import sys
import pysam
import numpy

#single test
#scripts/1.peak_shift.py data/6753_3-1-16.28116.IS614.bam tester tester2 tester3
#run with 
#ls -1 data/clinical | grep bam$ | xargs -n 1 -I foo sh -c "scripts/1.peak_shift.py data/clinical/foo output/peakshift_outputs/foo.corrected_depths.tsv output/peakshift_outputs/foo.cross_correlations.txt output/peakshift_outputs/foo.correctedreads.bed output/peakshift_outputs/foo.corrected.bam &"  
#ls -1 data/simulated_insertions | grep bam$ | xargs -n 1 -I foo sh -c "scripts/1.peak_shift.py data/clinical/foo output/peakshift_outputs/foo.corrected_depths.tsv output/peakshift_outputs/foo.cross_correlations.txt output/peakshift_outputs/foo.correctedreads.bed output/peakshift_outputs/foo.corrected.bam &"  


#then sort/index bam outputs
#ls -1 output/peakshift_outputs/*corrected.bam | xargs -n 1 -I foo sh -c "samtools sort -o foo.sorted.bam foo &"
#rm output/peakshift_outputs/*corrected.bam
#ls output/peakshift_outputs/*.bam.corrected.bam.sorted.bam | sed 's/.bam.corrected.bam.sorted.bam//g' | xargs -n 1 -I foo mv foo.bam.corrected.bam.sorted.bam foo.bam
#ls -1 output/peakshift_outputs/*.bam | xargs -n 1 -I foo sh -c "samtools index foo &"

def main():

	#IO
	if len(sys.argv) != 6:
		sys.exit('Usage: 1.peak_shift.py sample.bam corrected_depths.tsv cross_correlations_out.txt reads.bed shifted.bam')
		
	inputf = sys.argv[1]
	samfile_in = pysam.AlignmentFile(inputf, 'rb')		
	depths_outf = sys.argv[2]
	depths_out = open(depths_outf, 'w')
	cors_outf = sys.argv[3]
	cors_out = open(cors_outf, 'w')
	bed_outf = sys.argv[4]
	bed_out = open(bed_outf, 'w')
	bam_outf = sys.argv[5]
	bam_out = pysam.AlignmentFile(bam_outf, 'wb', template=samfile_in)
	
	
	#vars
	fwd_depths = []
	rev_depths = []
	chrom = ''
	start_idx = -1
	
    #extract coverages
	print("Counting stranded coverage")
	for pup in samfile_in.pileup(): #walk across one bp at a time
		fwd = [r for r in pup.pileups if not r.alignment.is_reverse]
		fwd_depths.append(len(fwd))
		rev = [r for r in pup.pileups if r.alignment.is_reverse]
		rev_depths.append(len(rev))
	
		if len(fwd) + len(rev) != pup.n: #sanity check
			sys.exit("something doesn't add up")
	
	print("Calculating cross-correlations")
	#check cross-correlations
	cors = []
	for offset in range(1,1000):
		offset_fwd = fwd_depths[0:(len(fwd_depths)-offset)]
		offset_rev = rev_depths[offset:len(rev_depths)]
		cor = numpy.correlate(offset_fwd, offset_rev)[0]
		cors_out.write(str(cor) + "\n")
		cors.append(cor)
	
	print("Reticulating splines")
	#bullshit
	print("Compiling modulo buffers")
	#more bullshit
	
	print("Generating output")
	
	print("Outputting coverages")
	#generate final output
	best_offset = numpy.argmax(cors)
	best_offset_fwd = fwd_depths[0:(len(fwd_depths)-best_offset)]
	best_offset_rev = rev_depths[best_offset:len(rev_depths)]
	
	chrom = samfile_in.get_reference_name(0)
	
	for idx in range(0, len(best_offset_fwd)):
		nextline = "\t".join([chrom, str(idx), str(best_offset_fwd[idx]+best_offset_rev[idx])]) + "\n"
		depths_out.write(nextline)
		
	
	print("Outputting shifted read BED and BAM files")
	for read in samfile_in.fetch(): #walk across one bp at a time
		if read.is_reverse:
			bed_out.write("\t".join([read.query_name, chrom, str(read.get_reference_positions()[0] - best_offset), str(max([0,read.get_reference_positions()[-1] - best_offset])), 'rev', str(-1*best_offset)]) + "\n")
			read.reference_start = max([0,read.get_reference_positions()[-1] - best_offset])
		else:
			bed_out.write("\t".join([read.query_name, chrom, str(read.get_reference_positions()[0] + best_offset), str(read.get_reference_positions()[-1] + best_offset), 'fwd', str(best_offset)]) + "\n")
			read.reference_start = read.reference_start + best_offset
		bam_out.write(read)
	samfile_in.close()
	cors_out.close()
	bam_out.close()
	
	
if __name__ == '__main__':
	main()