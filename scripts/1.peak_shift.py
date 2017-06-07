#!/usr/bin/env python

import sys
import pysam
import numpy

#single test
#scripts/1.peak_shift.py data/clinical/6753_3-1-16.28116.IS614.bam tester tester2 tester3 

#run with 
'''
ls -1 data/clinical | grep bam$ | xargs -n 1 -I foo sh -c "scripts/1.peak_shift.py data/clinical/foo  output/peakshift_outputs/foo.cross_correlations.txt output/peakshift_outputs/foo.offset.txt output/peakshift_outputs/foo.corrected.bam &"  
ls -1 data/simulated_insertions | grep bam$ | xargs -n 1 -I foo sh -c "scripts/1.peak_shift.py data/simulated_insertions/foo output/peakshift_outputs/foo.cross_correlations.txt output/peakshift_outputs/foo.offset.txt output/peakshift_outputs/foo.corrected.bam &"

ls -1 data/simulated_differential_insertions | grep bam$ | xargs -n 1 -I foo sh -c "scripts/1.peak_shift.py data/simulated_differential_insertions/foo  output/peakshift_outputs/foo.cross_correlations.txt output/peakshift_outputs/foo.offset.txt output/peakshift_outputs/foo.corrected.bam &"  
'''
#then sort/index bam outputs
#ls -1 output/peakshift_outputs/*corrected.bam | xargs -n 1 -I foo sh -c "samtools sort -o foo.sorted.bam foo &"
#rm output/peakshift_outputs/*corrected.bam
#ls output/peakshift_outputs/*.bam.corrected.bam.sorted.bam | sed 's/.bam.corrected.bam.sorted.bam//g' | xargs -n 1 -I foo mv foo.bam.corrected.bam.sorted.bam foo.bam
#ls -1 output/peakshift_outputs/*.bam | xargs -n 1 -I foo sh -c "samtools index foo &"

def main():
	if len(sys.argv) != 5:
		sys.exit('Usage: 1.peak_shift.py sample.bam cross_correlations_out.txt offset.txt shifted.bam')
		
	#open files named in parameters
	inputf = sys.argv[1]
	samfile_in = pysam.AlignmentFile(inputf, 'rb')		
	cors_outf = sys.argv[2]
	cors_out = open(cors_outf, 'w')
	offset_outf = sys.argv[3]
	offset_out = open(offset_outf, 'w')
	bam_outf = sys.argv[4]
	bam_out = pysam.AlignmentFile(bam_outf, 'wb', template=samfile_in)
	
	fwd_depths = []
	rev_depths = []
	
	#extract coverages by read orientation
	print("Counting stranded coverage")
	for pup in samfile_in.pileup(): #walk across positions with aligned reads
		fwd = [r for r in pup.pileups if not r.alignment.is_reverse]
		fwd_depths.append(len(fwd))
		rev = [r for r in pup.pileups if r.alignment.is_reverse]
		rev_depths.append(len(rev))
		
		if len(fwd) + len(rev) != pup.n: #sanity check
			sys.exit("something doesn't add up")

	#check cross-correlations
	print("Calculating cross-correlations")
	cors = numpy.correlate(fwd_depths, rev_depths, "same") #generate cross-correlations for one-sided overlaps of the two vectors of depths
	cors_out.write("\n".join([str(i) for i in cors])) #store the cross-correlations to a file
	#best offset is given by the degree of overlap between the 
	#two matrices at each cross-correlation calculation. 
	best_offset = len(cors)/2 - numpy.argmax(cors) 
	best_offset = best_offset/2 #divide by two to get the amount by which 
								#forward and reverse reads should each be shifted
	
	
	print("Outputting shifted read BAM files")
	for read in samfile_in.fetch(): #walk across one bp at a time
		if read.is_reverse:			
			read.reference_start = max([0,read.reference_start - best_offset])
		else:
			read.reference_start = read.reference_start + best_offset
			
		bam_out.write(read)
	
	#store the best offset to a file
	offset_out.write(str(best_offset))
	
	samfile_in.close()
	cors_out.close()
	offset_out.close()
	bam_out.close()
	
	
if __name__ == '__main__':
	main()