import sys
from snakemake import shell
from contextlib import redirect_stdout
from io import StringIO

def main(peak_files):

    merged_peaks = shell('cat %s | bedtools sort -i stdin | bedtools merge -i stdin' % ' '.join(peak_files),
                         iterable=True)


    for peak in merged_peaks:
        print(peak.split())



if __name__ == '__main__':
    peak_files = sys.argv[1:]

    main(peak_files)