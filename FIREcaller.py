"""Detect frequently interacting regions (FIREs) from Hi-C data. Ported from https://github.com/yycunc/FIREcaller."""

import argparse
import logging
import pandas as pd
import cooler
import statsmodels.api as sm
import numpy as np
import scipy.stats as st

logging.basicConfig(level=logging.DEBUG)

def count_cis_neighbors(mat, hic, bin_size):
    """For each region from mappability file return the number of cis-neighbors within given range.

    Parameters
    ----------
    mat : DataFrame
        Pandas DataFrame consisting of genomic regions
    hic : Cooler
        HiC experiment of size bin_size
    bin_size : int
        Bin size
    """
    chr_size = {hic.chromnames[i]:hic.chromsizes[i] for i in range(len(hic.chromnames))}

    cis_neighbours = list()
    i = 0
    for row in mat.itertuples(index=True):
        if row.chr in chr_size.keys() and row.end < chr_size[row.chr]:
            # calculate left cis-neighbors
            start = max(row.start - args.neighborhood_region, 0)
            end = max(row.end - bin_size, 0)
            cis_left = hic.matrix(balance=False).fetch(f"{row.chr}:{row.start}-{row.end}", f"{row.chr}:{start}-{end}")

            # calculate right cis-neighbors
            start = min(row.start + bin_size , chr_size[row.chr])
            end = min(row.end + args.neighborhood_region, chr_size[row.chr])
            cis_right = hic.matrix(balance=False).fetch(f"{row.chr}:{row.start}-{row.end}", f"{row.chr}:{start}-{end}")

            cis_neighbor_count = sum(cis_left[0]) + sum(cis_right[0])
        else:
            cis_neighbor_count = 0
        cis_neighbours.append(cis_neighbor_count)

        i+=1
        if i % 1000 == 0:
            logging.debug(f"{100*i/len(mat):.2f}%")
    return cis_neighbours

def remove_bad_regions(mat, bin_no, perc_threshold, avg_mappability_threshold):
    """Remove "bad" regions.

    Parameters:
    ----------
    mat : DataFrame
        Pandas DataFrame consisting of genomic regions
    bin_no : int
        number of adjacent bins considered as neighbors
    perc_threshold : float
        maximum ratio of "bad" neighbors allowed
    avg_mappability_threshold : float
        minimum mappability allowed
    """
    # temporary column names
    flag, bad_neig, perc = "flag", "bad_neig", "perc"

    # flag zeros
    mat[flag] = (mat["F"]==0) | (mat["GC"]==0) | (mat["M"]==0)

    # count "bad" (flagged) cis-neighbors
    mat[bad_neig] = 0
    for i in range(len(mat)):
        mat.at[i, bad_neig] = mat[flag].iloc[max(i-bin_no,0):max(i-1,0)].sum() + mat[flag].iloc[min(i+1,len(mat)):min(i+bin_no,len(mat))].sum()

    # calculate % of "bad" neighbors
    mat[perc] = mat[bad_neig] / (2*bin_no)

    # remove flag==1 and bad_neig <= perc_threshold && M > avg_mappability_threshold
    mat = mat[ (mat[flag] == False) & (mat[perc] <= perc_threshold) & (mat["M"] > avg_mappability_threshold)]

    # drop no-longer used columns
    mat = mat.drop(columns=[flag, bad_neig, perc])

    return mat

def hic_norm(mat, count_neig, fire):
    """Poisson normalization.

    Parameters:
    ----------
    mat : DataFrame
        Pandas DataFrame consisting of genomic regions
    count_neig : str
        DataFrame column name where neighbor count is located
    fire : str
        DataFrame column name where fire score should be stored
    """
    y = mat[count_neig]
    x = mat[["F", "GC", "M"]]
    x = sm.add_constant(x)

    glm = sm.GLM(y, x, family=sm.families.Poisson())
    res = glm.fit()        
    mat[fire] = mat[count_neig] / np.exp(res.params[0]+mat["F"]*res.params[1]+mat["GC"]*res.params[2]+mat["M"]*res.params[3])

def quantile_normalize(df_input, columns=None):
    """Perform quantile normalization on the selected columns from a DataFrame.

    Adopted from: https://github.com/ShawnLYU/Quantile_Normalize/blob/master/quantile_norm.py

    Parameters:
    ----------
    df_input : DataFrame
        DataFrame to normalize
    columns : list()
        List of column names to perform the normalization on
    """
    df = df_input.copy()
    if not columns:
        columns = df.columns
    #compute rank
    dic = {}
    for col in columns:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in columns:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df


def fire_caller(mat, fire, logpval):
    """Calculate log p-value.

    Parameters:
    ----------
    mat : DataFrame
        Pandas DataFrame consisting of genomic regions
    fire : str
        DataFrame column name where fire is located
    logpval : str
        DataFrame column name where log p-value should be stored
    """
    fire_mean = mat[fire].mean()
    fire_std = mat[fire].std()

    log_pvalue = - np.log( 1 - st.norm.cdf(mat[fire], loc=fire_mean, scale=fire_std) )
    mat[logpval] = log_pvalue

def calc_fires(mappability_filename, cooler_filenames, bin_size, neighborhood_region, perc_threshold=.25, avg_mappability_threshold=0.9):
    """Perform FIREcaller algorithm.

    Parameters:
    ----------
    mappability_filename : str
        Path to mappability file
    cooler_filenames : str
        List of paths to HiC experiment results in cooler format
    bin_size : int
        Bin size
    neighborhood_region : int
        Size of neighbor region
    perc_threshold : float
        maximum ratio of "bad" neighbors allowed
    avg_mappability_threshold : float
        minimum mappability allowed        
    """
    bin_size = bin_size
    bin_no = neighborhood_region // bin_size

    # read the mappability file
    mat = pd.read_csv(mappability_filename, delim_whitespace=True)
    required_cols = ["chr","start","end","F","GC","M"]
    if not all(col in mat.columns for col in required_cols):
        print(f"Error: Mappability file: {mappability_filename} does not contain all the required columns: {','.join(required_cols)}")
        exit()
    mp_bin_size = mat["end"].iloc[0] - mat["start"].iloc[0]
    if mp_bin_size != bin_size:
        print(f"Error: Mappability file: {mappability_filename} bin size is {mp_bin_size} while --bin_size={bin_size} was provided")
        exit()

#    mat = mat.iloc[:3000] # uncomment to limit the matrix size. useful while debugging

    for c, cooler_filename in enumerate(cooler_filenames):
        hic = cooler.Cooler(f"{cooler_filename}::resolutions/{bin_size}")
        # count cis neighbors
        mat[f"{c}_count_neig"] = count_cis_neighbors(mat, hic, bin_size)

    # remove "bad" regions
    mat = remove_bad_regions(mat, bin_no, perc_threshold, avg_mappability_threshold)

    # HiCNormCis 
    for c, _ in enumerate(cooler_filenames):
        hic_norm(mat, f"{c}_count_neig", f"{c}_fire")

    #quantile normalization
    if len(cooler_filenames)>1:
        fire_columns = [f"{c}_fire" for c, _ in enumerate(cooler_filenames)]
        mat = quantile_normalize(mat, fire_columns)

    ### fire caller
    for c, _ in enumerate(cooler_filenames):
        fire_caller(mat, f"{c}_fire", f"{c}_logpvalue")

    return mat

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mappability_filename", help="Mappability_file", required=True)
    parser.add_argument("--cooler_filenames", help="Cooler files (.mcool) to perform FIRE on", action='append', required=True)
    parser.add_argument("--output_filename", help="Output filename", required=True)
    parser.add_argument("--bin_size", help="Bin size", type=int, required=True)
    parser.add_argument("--neighborhood_region", help="The size of the cis-neighborhood region", type=int, default=200000)
    parser.add_argument("--perc_threshold", help="Maximum ratio of \"bad\" neighbors allowed", type=float, default=.25)
    parser.add_argument("--avg_mappability_threshold", help="Minimum average mappability allowed", type=float, default=0.9)
    args = parser.parse_args()

    mat = calc_fires(args.mappability_filename, args.cooler_filenames, args.bin_size, args.neighborhood_region, args.perc_threshold, args.avg_mappability_threshold)

    mat.to_csv(args.output_filename, sep=" ", index=False)