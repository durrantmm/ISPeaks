import sys
import pysam
from scipy.stats import poisson
import math

#Usage: scripts/sig_peaks.py output/peakshift_outputs/sim_bcaccae1.47678.random_seq1.bam.corrected test_peak_shift.tsv peakshift_length

window_size = 0
final_peaks = []
def main():
	if len(sys.argv) != 4:
		print("error incorrect arg length")
		return
	global window_size
	input_file = sys.argv[1]
	output_file = sys.argv[2]
	window_size = 2 * int(sys.argv[3])

	samfile = pysam.AlignmentFile(input_file, 'rb')
	for contig in samfile.header['SQ']: #go through each contig to process peaks within it
		contig_name = contig['SN'] #get the name of the contig
		contig_length = contig['LN'] #get the length of the contig

		total_region = contig_name.strip()
		total_count = samfile.count(region=total_region)
		lambda_bg = (float(total_count) / contig_length) * window_size #lambda background calculated per window unit

		window_arr = [] #array of windows

		for pup in samfile.pileup(region=contig_name): #do calculations for each window (sliding base by base over regions where we have read coverage)
			start_pos = pup.pos
			end_pos = pup.pos + window_size
			depth = float(samfile.count(contig_name, start_pos, end_pos))
			p_value = (1 - poisson.cdf(depth, mu=lambda_bg)) #calculate p_value from the poisson
			window_obj = {"start":start_pos, "depth": depth, "p_value": p_value}
			print(contig_name, window_obj)
			window_arr.append(window_obj)

		p_thresh = 0.00005 #set the p-threshold, IMPORTANT
		active_peak = False
		current_peak = {}
		peaks = []
		for window in window_arr: #merge the overlapping significant windows below
			if window["p_value"] <= p_thresh: #found a significant window
				if not active_peak:
					current_peak = window
					current_peak["extension"] = window["start"]
					current_peak["contig"] = contig_name
					active_peak = True
				else: #potentially overlapping the current peak
					curr_end = current_peak["start"] + window_size
					if window["start"] >= current_peak["start"] and window["start"] <= curr_end:
						current_peak["p_value"] = min(current_peak["p_value"], window["p_value"]) #get lower p-value
						current_peak["extension"] = window["start"] #merge the peaks
					else:
						peaks.append(current_peak) #not overlapping anymore -> current peak is finished
						current_peak = {}
						active_peak = False
			else:
				if active_peak:
					peaks.append(current_peak)
					current_peak = {}
					active_peak = False
		process_final_peaks(peaks, samfile) #do final peak processing

	output_final(output_file , final_peaks) #print final bed output



def process_final_peaks(peaks, samfile):
	for peak in peaks:
		cont_name = peak["contig"] #set peak name
		start_pos = peak["start"] #set peak start position
		end_pos = peak["extension"] + window_size #set peak end position
		peak["end"] = end_pos #set peak end position
		peak["peak_name"] = cont_name + "_" + str(start_pos) #name the peak
		region_name = cont_name.strip() + ":" + str(start_pos) + ":" + str(end_pos) #define the region of the peak
		peak["count"] = samfile.count(region=region_name) #count reads within the region
		if peak["p_value"] == 0: #correct so the log does not crash the program
			peak["p_value"] = 0.000000001 #small unimportant value
		peak["log_p"] = (- math.log10(peak["p_value"])) #set the -log10(p_value) of the peak
		pileup_depths = [len(pup.pileups) for pup in samfile.pileup(region=region_name)] #search for peak within region
		if len(pileup_depths)== 0:
			max_depth = 0
		else:
			max_depth = max(pileup_depths)
		peak["peak"] = start_pos + max_depth #set the peak peak
		final_peaks.append(peak) #add to the final_peaks global array



def output_final(f_name, peaks):
	print("Outputting data")
	bed_out = open(f_name, "w")
	for p in peaks: #write a new line for each of the peaks
		bed_out.write("\t".join([p["contig"], str(p["start"]), str(p["end"]), p["peak_name"], "0", "+", str(p["count"]), str(p["log_p"]), "-1", str(p["peak"])]) + "\n")

	bed_out.close()


if __name__ == '__main__':
	main()
