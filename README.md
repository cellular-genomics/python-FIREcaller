# Detect Frequently Interacting REgions (FIREs) in Python

[![PyPI](https://img.shields.io/pypi/v/FIREcaller.svg)](https://pypi.python.org/pypi/FIREcaller)

The project is a port of [the R package for detecting frequently interacting regions (FIREs) from Hi-C data](https://github.com/yycunc/FIREcaller) to Python. FIRE is described in [A Compendium of Chromatin Contact Maps Reveal Spatially Active Regions in the Human Genome](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5478386/) paper.

## Command line usage

Install the package:
```
python3 -m pip install FIREcaller
```

Download [some HiC experiment results](https://data.4dnucleome.org/browse/?experiments_in_set.biosample.biosource.individual.organism.name=human&experiments_in_set.experiment_type.display_title=in+situ+Hi-C&experimentset_type=replicate&type=ExperimentSetReplicate) from the 4D Nucleome database. Choose Contact Matrix (.mcool) as a file type.

Download the mappability file from [Yunjiang's website](http://enhancer.sdsc.edu/yunjiang/resources/genomic_features/). Ensure the genomic assembly and the bin size matches your HiC files. Add header line `chr start end F GC M` if missing.

Perform FIRE calling:
```
FIREcaller \
  --cooler_filenames 4DNFIT5YVTLO.mcool \
  --cooler_filenames 4DNFIJTOIGOI.mcool \
  --mappability_filename F_GC_M_Hind3_10Kb_el.GRCh38.txt \
  --bin_size 10000 \
  --output_filename fires.csv
```

The output file will consist of the genomic regions and their corresponding FIRE scores and log p-values for each HiC file:
```
chr start end F GC M 0_count_neig 1_count_neig 0_fire 1_fire 0_logpvalue 1_logpvalue
chr1 1970000 1980000 2000 0.5175 1.0000 5365 2376 0.9515 0.8702 0.5861 0.4317
chr1 2020000 2030000 3000 0.6287 0.9907 4305 2005 0.5806 0.5831 0.1128 0.1144
chr1 2060000 2070000 4000 0.4770 0.9210 4029 2171 0.6880 0.8678 0.1954 0.4277
(...)
```

`{n}_fire` column stores the FIIRE score for the n-th cooler file provided as `--cooler_filenames` argument

`{n}_logpvalue` column stores the log p-value for the n-th cooler file provided as `--cooler_filenames` argument

## Python usage

Use the `FIREcaller.calc_fires()` function to perform FIRE calling from your Python program:

`calc_fires(mappability_filename, cooler_filenames, bin_size, neighborhood_region, perc_threshold=.25, avg_mappability_threshold=0.9)`

`mappability_filename : str` - Path to mappability file

`cooler_filenames : str` -  List of paths to HiC experiment results in cooler format

`bin_size : int` - Bin size

`neighborhood_region : int` - Size of neighbor region

`perc_threshold : float` - Maximum ratio of "bad" neighbors allowed

`avg_mappability_threshold : float` - Minimum mappability allowed        

The function returns the Pandas DataFrame matrix.

## Verification

In order to compare the results derived from this project with the ones obtained from the [FIREcaller R package](https://github.com/yycunc/FIREcaller/) please follow the steps described in the [`README.md`](https://github.com/yycunc/FIREcaller/blob/master/README.md) file to produce the `FIRE_ANALYSIS_40KB.txt` file.

Now download and uncompress the [Hippo_Hi-C_inputs_chr1_22.tar.gz](https://yunliweb.its.unc.edu/FIREcaller/example/HiC_input_for_FIREcaller/Hippo_Hi-C_inputs_chr1_22.tar.gz) and [F_GC_M_HindIII_40KB_hg19.txt.gz](https://yunliweb.its.unc.edu/FIREcaller/example/HiC_input_for_FIREcaller/F_GC_M_HindIII_40KB_hg19.txt.gz) files.

Convert the HiC files to get the `hippo.mcool` cooler file using the `scripts/build_hippo_cool.py` python script:

```
python scripts/build_hippo_cool.py
```

Run the FIRE calling:
```
FIREcaller \
  --cooler_filename hippo.mcool \
  --mappability_filename F_GC_M_HindIII_40KB_hg19.txt \
  --bin_size 40000 \
  --output_filename hippo_fires.csv
```

Compare the `FIRE_ANALYSIS_40KB.txt` and `fires.csv` files.

## Citation

Crowley, C., Yang, Y.*, Qiu, Y., Hu, B., Abnousi, A., Lipiński, J., Plewczynski, D., Wu, D., Won, H., Ren, B., Hu, M.*, Li, Y*. (2021) FIREcaller: Detecting Frequently Interacting Regions from Hi-C Data. Computational and Structural Biotechnology Journal, 19: 355–362.

## TODO:

- [ ] Add option to remove the MHC regions
- [ ] Calculate the Super FIREs
- [x] Convert the script into the Python package

Let me know if you need any of the above and I'll be happy to implement it. I also accept PRs and comments. Enjoy.
