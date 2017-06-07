import sys


def main(peaks_path, sites_path):

    peaks = load_peaks(peaks_path)
    sites = load_sites(sites_path)

    tp = count_true_positives(peaks, sites)
    fp = count_false_positives(peaks, sites)
    fn = count_false_negatives(peaks, sites)


    print("True Positives: %d" % tp)
    print("False Positives: %d" % fp)
    print("False Negatives: %d" % fn)

def count_true_positives(peaks, sites):

    count = 0
    for peak in peaks:
        p_contig, start, end = peak
        for site in sites:
            s_contig, loc = site
            if p_contig == s_contig:
                if start <= loc <= end:
                    count += 1
    return count

def count_false_positives(peaks, sites):
    count = 0
    for peak in peaks:
        p_contig, start, end = peak
        true_pos = False
        for site in sites:
            s_contig, loc = site
            if p_contig == s_contig:
                if start <= loc <= end:
                    true_pos = True
        if not true_pos:
            count += 1
    return count

def count_false_negatives(peaks, sites):
    count = 0
    for site in sites:
        s_contig, loc = site
        true_pos = False
        for peak in peaks:
            p_contig, start, end = peak
            if start <= loc <= end:
                true_pos = True
        if not true_pos:
            count += 1

    return count


def load_peaks(peaks_path):
    peaks = []
    with open(peaks_path) as infile:
        for line in infile:
            line = line.strip().split()
            peak = line[0], int(line[1]), int(line[2])
            peaks.append(peak)
    return peaks


def load_sites(sites_path):
    sites = []
    with open(sites_path) as infile:
        for line in infile:
            line = line.strip().split()
            site = line[0].split(':')
            site = site[0], int(site[1])
            sites.append(site)
    return sites


if __name__ == '__main__':

    peaks = sys.argv[1]
    sites = sys.argv[2]

    main(peaks, sites)