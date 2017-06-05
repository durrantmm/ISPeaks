import sys
from scipy.stats import poisson
import numpy as np


#Usage scripts/2.significant_peaks.py output/peakshift_outputs/sim_bcaccae4.47678.random_seq1.bam.corrected_depths.tsv

def main():
	read_f = sys.argv[1] #input corrected_depths file
	out_f = sys.argv[2] #peak output file (not used yet..)

	max_length = {}
	data_obj = {}

	#Get the maximum position size to create empty arrays
	with open(read_f, "r") as reads:
		for line in reads:
			chrom, pos, count = line.strip().split('\t')
			if chrom not in max_length:
				max_length[chrom] = 0
			else:
				if pos > max_length[chrom]:
					max_length[chrom] = int(pos)
	print(max_length) #maximum bp position of each contig

	for chrom in max_length: #instantiate a data object for each chromosome/contig
		data_obj[chrom] = {"reads": [], "p_vals": [], "lams": [], "total_count": 0, "total_length":max_length[chrom], "reads_with_counts":0}
		data_obj[chrom]["reads"] = [0 for x in range(max_length[chrom] + 2)]

	# Update each position with the reads from the corrected_depths files
	with open(read_f, "r") as reads:
		for line in reads:
			chrom, pos, count = line.strip().split('\t')
			if int(count) > 0:
				data_obj[chrom]["reads_with_counts"] += 1
			data_obj[chrom]["reads"].insert(int(pos), int(count))
			data_obj[chrom]["total_count"] += int(count)

	# Calculate p-values for each base
	for chrom in data_obj:
		lambda_loc = data_obj[chrom]["total_count"] / data_obj[chrom]["reads_with_counts"]
		#single lambda_loc value used to calculate all p-values
		print("count: ", data_obj[chrom]["total_count"] )
		print("length: ", data_obj[chrom]["total_length"] )
		print("lambda", lambda_loc)
		#Calculate each p-value for each position
        data_obj[chrom]["p_vals"] = [(1 - poisson.cdf(r, mu=lambda_loc)) for r in data_obj[chrom]["reads"]]

    #Detect peaks (contiguous significant bases)
	for chrom in data_obj:
		raw_peaks = [] #array of peak objects
		peak_active = False
		missed_bases = 0
		current_peak = {}
		p_vals = data_obj[chrom]["p_vals"]
		for pos in range(len(p_vals)):
			if p_vals[pos] < 0.0005: #p-value threshold (can update with q-value)
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
					peak_active = False #not currently adding to a peak
		print(raw_peaks) #unprocessed peaks
		print("total raw peaks: ", len(raw_peaks)) 
		processed_peaks = []
		for peak in raw_peaks:
			length = peak['end'] - peak['start']
			if length >= 10: #only process peaks larger than a given length (10)
				peak["length"] = length
				depth = peak['total_depth']
				peak["coverage"] = depth / length
				processed_peaks.append(peak)


		print("total processed peaks: ", len(processed_peaks)) #processed peak length
		print("processed peaks: ", processed_peaks) #processed peaks



# Update the parameter peak object (new endpoint and add to total depth)
def update_peak_obj(peak, pos, reads):
	peak["end"] = pos 
	peak["total_depth"] += reads
	#do other stuff for the peak here
	return peak 


# Return a newly created "peak" object
def create_peak_obj(start_pos, reads):
	peak = {"start": start_pos, "end": start_pos, "total_depth": reads}
	return peak 
	#add other stats about peak...




if __name__ == '__main__':
    main()


