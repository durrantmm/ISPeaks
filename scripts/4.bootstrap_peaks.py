import sys
from collections import defaultdict
import pysam
from snakemake import shell
import numpy as np

def main(peaks_path, bam_path, bootstraps=10):

    peaks = load_peaks(peaks_path)
    bam = pysam.AlignmentFile(bam_path, 'rb')

    total_reads = sum([1 for read in bam])

    read_counts = dict()

    emp_read_counts = empirical_read_counts(peaks, bam)
    read_counts['EMPIRICAL'] = emp_read_counts

    for b in range(bootstraps):
        read_counts['BOOTSTRAP%d' % (b+1)] = bootstrap_read_counts(peaks, bam_path, total_reads)

    print('\t'.join(['Sampling', 'Contig', 'PeakStart', 'PeakEnd', 'ReadCount', 'ReadDepth']))
    for iter in read_counts:
        for contig in read_counts[iter]:
            for peak in read_counts[iter][contig]:
                print('\t'.join([iter, contig, '\t'.join(map(str, peak))]))


def bootstrap_read_counts(peaks, bam_path, total_reads):

    indices = sorted(list(np.random.choice(range(total_reads), total_reads, replace=True)))

    i = 0

    bam = pysam.AlignmentFile(bam_path, 'rb')
    outfile = pysam.AlignmentFile("boot.tmp.sam", "w", template=bam)

    for read in bam:

        try:
            while i == indices[0]:
                outfile.write(read)
                indices = indices[1:]
        except:
            break
        i += 1

    outfile.close()

    shell('samtools view -b boot.tmp.sam -o boot.tmp.bam')
    shell('samtools sort -o boot.tmp.sorted.bam boot.tmp.bam')
    shell('samtools index boot.tmp.sorted.bam')

    boot_read_counts = empirical_read_counts(peaks, pysam.AlignmentFile('boot.tmp.sorted.bam', 'r'))
    shell('rm boot.tmp*')

    return boot_read_counts


def empirical_read_counts(peaks, bam):
    out_data = defaultdict(list)

    for contig in peaks:
        for peak in peaks[contig]:
            count = bam.count(contig, peak[0], peak[1])
            pileup_depths = [len(pup.pileups) for pup in bam.pileup(contig, peak[0], peak[1])]

            if len(pileup_depths)== 0:
                max_depth = 0
            else:
                max_depth = max(pileup_depths)

            out_data[contig].append(peak + [count, max_depth])

    return out_data


def load_peaks(peaks_path):
    peaks = defaultdict(list)
    with open(peaks_path) as infile:
        for line in infile:
            line = line.strip().split()
            contig, start, end = line[0], int(line[1]), int(line[2])
            peaks[contig].append([start, end])

    return dict(peaks)


if __name__ == '__main__':

    merged_peaks_path = sys.argv[1]
    bam = sys.argv[2]

    main(merged_peaks_path, bam)