import sys
from snakemake import shell
from contextlib import redirect_stdout
from io import StringIO

usage = "python scripts/3.merge_peaks.py peaks.in.1.tsv peaks.in.2.tsv peaks.in.3.tsv ... > peaks.out.tsv"

def main(peak_files):

    # Makes a call to the shell environment to use bedtools to carry out the required operation
    merged_peaks = shell('cat %s | bedtools sort -i stdin | bedtools merge -i stdin' % ' '.join(peak_files),
                         iterable=True)


    for peak in merged_peaks:
        print(peak)



if __name__ == '__main__':
    # A list of peak files to merge, in bed format.
    if len(sys.argv) == 1:
        print(usage)
        sys.exit()

    peak_files = sys.argv[1:]

    main(peak_files)