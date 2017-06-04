# USAGE: python scripts/2.significant_peaks_mdurrant.py output/peakshift_outputs/sim_bcaccae4.47678.random_seq1.bam.correctedreads.bed data/47678.fasta


import sys
from Bio import SeqIO
import numpy as np
from scipy.stats import poisson
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import pandas as pd

def main(bedpath, genome):

    # Processing genome
    genome_arr = process_genome(genome)

    # Adding read counts
    genome_arr = add_read_counts(genome_arr, bedpath)

    # Calling peaks
    peaks = call_peaks(genome_arr)

    for peak in peaks:
        print('\t'.join(map(str, peak)))

def call_peaks(genome, unit_length = 200, small_length = 1000, medium_length = 5000, large_length = 10000):

    peaks_out = []
    for contig in genome:

        total_reads = sum(genome[contig])
        contig_length = len(genome[contig])

        if total_reads == 0:
            continue

        window_counts = window_read_counts(genome[contig], unit_length)
        window_sum = window_counts.sum(axis=1)

        small_bin_counts = calculate_bins(window_counts, small_length, unit_length)
        medium_bin_counts = calculate_bins(window_counts, medium_length, unit_length)
        large_bin_counts = calculate_bins(window_counts, large_length, unit_length)

        local_bin_sums = np.hstack((small_bin_counts.sum(axis=1),
                                    medium_bin_counts.sum(axis=1),
                                    large_bin_counts.sum(axis=1)))

        local_lambdas = (local_bin_sums / np.array([small_length, medium_length, large_length])) * unit_length
        lambda_bg = np.ones(window_sum.shape) * (total_reads / contig_length) * unit_length

        all_lambdas = np.hstack((lambda_bg, local_lambdas))
        max_lambdas = np.amax(all_lambdas, axis=1)

        p_vals = 1 - poisson.cdf(window_sum.astype(int), mu=max_lambdas.astype(float))
        p_vals = np.transpose(p_vals.astype(np.longdouble))[0]

        qvalue = importr('qvalue')
        q_vals = np.array(qvalue.qvalue(FloatVector(p_vals))[2])
        q_vals = np.hstack((np.transpose(np.matrix(list(range(1, len(q_vals) + 1)))),
                            np.transpose(np.matrix(q_vals))))
        qv_df = pd.DataFrame(q_vals)
        qv_df.columns = ['Position', 'qvalue']

        peak_indices = np.array(qv_df.query('qvalue < 0.01')['Position'].tolist()).astype(int)

        peaks = indices_to_peaks(peak_indices)
        peaks = correct_peaks(peaks, unit_length)

        for peak in peaks:
            peaks_out.append([contig, peak[0], peak[1]])

    return peaks_out


def indices_to_peaks(peak_indices):

    peaks = []
    new_peak = [0, 0]
    for i, e in enumerate(peak_indices):
        if new_peak == [0, 0]:
            new_peak[0] = e
        elif e - peak_indices[i - 1] > 1:
            new_peak[1] = peak_indices[i - 1]
            peaks.append(new_peak)
            new_peak = [e, 0]
    else:
        new_peak[1] = e
        peaks.append(new_peak)

    return peaks


def correct_peaks(peaks, unit_length):
    peaks_out = []
    for peak in peaks:
        peaks_out.append([peak[0] * unit_length, peak[1] * unit_length + unit_length])
    return peaks_out


def calculate_bins(window_counts, bin_size, unit_length):

    window_arr = np.array(np.reshape(window_counts,
                            (1, window_counts.shape[0] * window_counts.shape[1])))[0]

    window_ranges_left = range(0, window_counts.shape[0] * unit_length, unit_length)
    window_ranges_right = range(unit_length, window_counts.shape[0] * unit_length + unit_length, unit_length)

    windows = np.hstack((np.reshape(window_ranges_left, (len(window_ranges_left), 1)),
                         np.reshape(window_ranges_right, (len(window_ranges_right), 1))))

    bins = []

    for i, w in enumerate(windows):
        bin_range = [int(w[0] - bin_size/2), int(w[1] + bin_size/2)]

        if bin_range[0] < 0 and bin_range[1] > len(window_arr) - 1:
            bin = window_arr[0: len(window_arr) - 1]
            bin = np.array([0] * abs(bin_range[0]) + bin + [0] * (bin_range[1] - len(window_arr) + 1))

        elif bin_range[0] < 0:
            bin = window_arr[0: bin_range[1]]
            bin = np.hstack((np.array([0] * abs(bin_range[0])), bin))

        elif bin_range[1] > len(window_arr) - 1:
            bin = window_arr[bin_range[0]: len(window_arr)]
            bin = np.hstack((bin, np.array([0] * (bin_range[1] - len(window_arr)))))

        else:
            bin = window_arr[bin_range[0]: bin_range[1]]

        bin = np.hstack( (bin[: int(bin_size/2)], bin[int(bin_size/2) + unit_length:]))

        bins.append(bin)

    return np.matrix(bins)


def window_read_counts(arr, window_size):
    n_windows = int(len(arr) / window_size)

    if n_windows * window_size != len(arr):
        new_len = n_windows * window_size + window_size
        extension = np.zeros(new_len - len(arr))
        arr = np.hstack((arr, extension))

    arr = np.reshape(arr, (int(len(arr) / window_size), window_size))

    return np.matrix(arr)


def add_read_counts(genome, bedpath):

    with open(bedpath) as infile:

        for line in infile:

            line = line.strip().split()
            contig, start, end, strand = line[1], int(line[2]), int(line[3]), line[4]

            con_length = len(genome[contig])

            if strand == 'rev' and start < 0:
                continue
            elif strand == 'fwd' and end > con_length:
                continue
            elif strand == 'rev':
                genome[contig][start] += 1
            elif strand == 'fwd':
                genome[contig][end] += 1

    return genome




def process_genome(genome):
    genome_out = dict()
    for contig in SeqIO.parse(genome, "fasta"):
        name = contig.name
        length = len(contig.seq)
        counts = np.zeros(length, dtype=np.int)
        genome_out[name] = counts
    return genome_out

if __name__ == '__main__':

    inbed = sys.argv[1]
    genome = sys.argv[2]

    main(inbed, genome)