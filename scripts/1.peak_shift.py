#!/usr/bin/env python

import sys
import pysam
import numpy

#run with ls -1 data/| grep bam$ | xargs -n 1 -I foo sh -c "scripts/1.peak_shift.py data/foo output/foo.corrected_depths.tsv output/foo.cross_correlations.txt &"  

def main():

	#IO
	if len(sys.argv) != 4:
		sys.exit('Usage: 1.peak_shift.py sample.bam corrected_depths.tsv cross_correlations_out.txt')
	inputf = sys.argv[1]
	samfile = pysam.AlignmentFile(inputf, 'rb')		
	depths_outf = sys.argv[2]
	depths_out = open(depths_outf, 'w')
	cors_outf = sys.argv[3]
	cors_out = open(cors_outf, 'w')

	#vars
	fwd_depths = []
	rev_depths = []
	chrom = ''
	start_idx = -1
	
    #extract coverages
    
    #SHIFT REVERSE MAPPING POSITIONS???
    
    
	for pup in samfile.pileup(): #walk across one bp at a time
		fwd = [r for r in pup.pileups if not r.alignment.is_reverse]
		fwd_depths.append(len(fwd))
		rev = [r for r in pup.pileups if r.alignment.is_reverse]
		rev_depths.append(len(rev))
	
		if len(fwd) + len(rev) != pup.n: #sanity check
			sys.exit("something doesn't add up")

	#check cross-correlations
	cors = []
	for offset in range(1,1000):
		offset_fwd = fwd_depths[0:(len(fwd_depths)-offset)]
		offset_rev = rev_depths[offset:len(rev_depths)]
		cor = numpy.correlate(offset_fwd, offset_rev)[0]
		cors_out.write(str(cor) + "\n")
		cors.append(cor)

	#generate final output
	best_offset = numpy.argmax(cors)
	best_offset_fwd = fwd_depths[0:(len(fwd_depths)-best_offset)]
	best_offset_rev = rev_depths[best_offset:len(rev_depths)]
	
	chrom = samfile.getrname()
	
	for idx in range(0, len(best_offset_fwd)):
		nextline = "\t".join(['chrom', str(idx), str(best_offset_fwd[idx]+best_offset_rev[idx])]) + "\n"
		depths_out.write(nextline)
		

	samfile.close()
	cors_out.close()

if __name__ == '__main__':
	main()