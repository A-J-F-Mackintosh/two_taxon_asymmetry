# two taxon asymmetry
This repository contains Mathematica and python notebooks associated with Mackintosh and Setter (2024), as well as a python script to estimate A<sub>m</sub> from a VCF file.

The Mathematica notebook shows how the generating function of genealogical branch lengths (Lohse et al. 2011) can be used to obtain (both marginal and conditional) expectations for the length of external branches, and therefore A<sub>m</sub>.

The python notebooks use coalescent simulations with msprime to estimate A<sub>m</sub> for arbitrary demographic histories, block sizes and recombination rates. Notebook 1 lays out the general approach, while notebooks 2-4 simulate specific demographic histories.

The script `estimate_Am.py` can be used to estimate A<sub>m</sub> from a (filtered) VCF file.

```
[Options]
    -v, --vcf <STR>                             VCF file
    -a, --a_samples <STR>                       A comma delimited list of samples from population A, e.g. sample_1,sample_2,sample_3
    -b, --b_samples <STR>                       A comma delimited list of samples from population B
    -j, --jackknife_blocksize <INT>             Number of SNPs in each block (aim for >20 blocks) [default: 100_000]
    -r, --rollover <STR>                        Allow jackknife blocks to roll over between sequences [default: False]
    -s, --span <INT>                            Look for het' mutations within this many bases of a hetAB
    -h, --help                                  Show this message
```

An example command would be:

`estimate_Am.py -v heliconius_20220202.vcf.gz -a ros.CJ2071.m -b chi.CJ565.m -j 100_000 -s 8 > heliconius_Am_8.out`

Note that value for `-s` is equivalent to block_size = (2 * value) + 1. For example, searching 4 bases either side of a hetAB mutation is equivalent to a block size of 9 bases. It is always recommended to calculate A<sub>m</sub> across a range of blocksize.

The script requires the modules `docopt`, `scikit-allel` and `numpy` which can all be installed through conda.


