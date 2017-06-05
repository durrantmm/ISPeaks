import sys
from scipy.stats import poisson
import numpy as np


#Usage scripts/2.significant_peaks.py output/peakshift_outputs/sim_bcaccae4.47678.random_seq1.bam.corrected_depths.tsv

def main():
	read_f = sys.argv[1]
	out_f = sys.argv[2]

	max_length = {}
	data_obj = {}

	with open(read_f, "r") as reads:
		for line in reads:
			chrom, pos, count = line.strip().split('\t')
			if chrom not in max_length:
				max_length[chrom] = 0
			else:
				if pos > max_length[chrom]:
					max_length[chrom] = int(pos)
	print(max_length)

	#create new objs...
	for chrom in max_length:
		data_obj[chrom] = {"reads": [], "p_vals": [], "lams": [], "total_count": 0, "total_length":max_length[chrom], "reads_with_counts":0}
		data_obj[chrom]["reads"] = [0 for x in range(max_length[chrom] + 2)]

	with open(read_f, "r") as reads:
		for line in reads:
			chrom, pos, count = line.strip().split('\t')
			if int(count) > 0:
				data_obj[chrom]["reads_with_counts"] += 1
			data_obj[chrom]["reads"].insert(int(pos), int(count))
			data_obj[chrom]["total_count"] += int(count)

	
	for chrom in data_obj:
		lambda_loc = data_obj[chrom]["total_count"] / data_obj[chrom]["reads_with_counts"]
		print("count: ", data_obj[chrom]["total_count"] )
		print("length: ", data_obj[chrom]["total_length"] )
		print("lambda", lambda_loc)

        data_obj[chrom]["p_vals"] = [(1 - poisson.cdf(r, mu=lambda_loc)) for r in data_obj[chrom]["reads"]]

    #peak calling...
	for chrom in data_obj:
		raw_peaks = [] #array of objects
		peak_active = False
		missed_bases = 0
		current_peak = {}
		p_vals = data_obj[chrom]["p_vals"]
		for pos in range(len(p_vals)):
			if p_vals[pos] < 0.0005: #significant base pair...
				missed_bases = 0
				reads = data_obj[chrom]["reads"][pos]
				if peak_active: #add this base to current peak
					current_peak = update_peak_obj(current_peak, pos, reads)
				else: #create a new peak
					current_peak = create_peak_obj(pos, reads)
					peak_active = True 
			else: #not a significant base pair
				missed_bases += 1
				if missed_bases == 10 and peak_active: #peak has ended
					raw_peaks.append(current_peak)
					peak_active = False
		print(raw_peaks)
		print("total raw peaks: ", len(raw_peaks))
		processed_peaks = []
		for peak in raw_peaks:
			length = peak['end'] - peak['start']
			if length >= 10:
				peak["length"] = length
				depth = peak['total_depth']
				peak["coverage"] = depth / length
				processed_peaks.append(peak)


		print("total processed peaks: ", len(processed_peaks))
		print("processed peaks: ", processed_peaks)






def update_peak_obj(peak, pos, reads):
	peak["end"] = pos 
	peak["total_depth"] += reads
	#do other stuff for the peak here
	return peak 



def create_peak_obj(start_pos, reads):
	peak = {"start": start_pos, "end": start_pos, "total_depth": reads}
	return peak 
	#add other stats about peak...




if __name__ == '__main__':
    main()


