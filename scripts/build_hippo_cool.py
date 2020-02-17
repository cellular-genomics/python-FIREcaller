import csv
import cooler
from pandas import DataFrame

def read_bins(bin_size):
    bins = {"chrom":list(), "start":list(), "end":list()}
    for ch_no in range(1, 23):
        ch = f"chr{ch_no}"
        start = 0
        with open(f"Hippo_{ch}") as hic_file:
            for line in hic_file:
                if line.strip():
                    bins["chrom"].append(ch)
                    bins["start"].append(start)
                    bins["end"].append(start+bin_size)
                    start += bin_size

    return DataFrame(data = bins, copy=True)

def pixel_iter():
    chr_bin = 0
    for ch_no in range(1, 23):
        ch = f"chr{ch_no}"
        print(ch, chr_bin)
        with open(f"Hippo_{ch}") as hic_file:
            counts = {"bin1_id":list(), "bin2_id":list(), "count":list()}
            i = 0
            for line in hic_file:
                if line.strip():
                    row = line.split("\t")
                    for j, cnt in enumerate(row):
                        if j>=i and cnt.strip() != "":
                            counts["bin1_id"].append(chr_bin+i)
                            counts["bin2_id"].append(chr_bin+j)
                            counts["count"].append(int(cnt))
                    i += 1
            yield DataFrame(data = counts, copy=True)
            chr_bin += i

bin_size = 40000

bins = read_bins(bin_size)
pixels = pixel_iter()

cooler.create_cooler(f"hippo.mcool::resolutions/{bin_size}", bins, pixels, ordered=True)