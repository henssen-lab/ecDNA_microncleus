#takes bam as input and returns log2 per bin mean cov
import pysam
import numpy as np
import pandas as pd
import argparse

#argparse
parser = argparse.ArgumentParser()
parser.add_argument('--file', type=str, required=True)
parser.add_argument('--binsize', type=int, required=True)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()

bam =args.file
bin_size =args.binsize
outpath =args.output


#takes chromosome length and create bins of chosen length, last bin contains the remaining bases,
# if readlength is not a multiple of binsize; 0 indexed as always
def make_bins(chrom_length, binsize):
    bins = []
    for i in range(-1, chrom_length, binsize):  # start at -1 to adhere to 0 indexing
        bins.append([i + 1, i + binsize])
    bins.append([i + 1, chrom_length - 1])
    return(bins)


#counts per bin coverages and sums them together, divide by binsize to make a mean per base coverage,
#log2 transformed, returns list of list with chr,start,end,number of bases covered, log2meanperbase
def calc_log2_per_bin_mean_cov(bamfile,bin_size,reference):
    log2_mean_per_base_cov_list = []
    for chr in list(reference.keys()):
        for bin in make_bins(reference[chr],bin_size):
            n_bases_covered = sum([x.nsegments for x in bamfile.pileup(chr, start=bin[0],end=bin[1],max_depth = 50000)])
            log2_mean_per_base_cov = np.log2( n_bases_covered / (bin[1]-bin[0]))
            log2_mean_per_base_cov_list.append([chr,bin[0],bin[1],n_bases_covered,log2_mean_per_base_cov])
    return(log2_mean_per_base_cov_list)

def main():
    bamfile = pysam.AlignmentFile(bam, "rb")
    reference = dict(zip(bamfile.references, bamfile.lengths))
    log2_mean_cov = calc_log2_per_bin_mean_cov(bamfile, bin_size, reference)
    df = pd.DataFrame(log2_mean_cov, columns=['chr', 'start', 'end', 'n_bases_covered', 'log2_mean_cov'])
    with open(outpath, 'w', newline='') as file:
        df.to_csv(file, header=True, index=False, sep="\t")

if __name__ == '__main__':
    main()