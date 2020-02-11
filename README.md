#  Dsuite
Preprint now on bioRxiv:  
Dsuite - fast D-statistics and related admixture evidence from VCF files  
Milan Malinsky  
bioRxiv 634477; doi: https://doi.org/10.1101/634477  
Data used in the manuscript:
[VCF file](http://cichlid.gurdon.cam.ac.uk/Malinsky_et_al_2018_LakeMalawiCichlids_scaffold_0.vcf.gz); [SETS.txt file](http://cichlid.gurdon.cam.ac.uk/sets.txt)

There is also a very nice and detailed [tutorial](https://github.com/mmatschiner/tutorials/tree/master/analysis_of_introgression_with_snp_data) by [@mmatschiner](https://github.com/mmatschiner).

## Compilation

To compile you must have a reasonably recent GCC (>=4.9.0) or clang compiler (on mac OS this comes with Command Line Tools) and the zlib compression library (https://www.zlib.net). Both should already be present on most systems. 


```console
$ git clone https://github.com/millanek/Dsuite.git
$ cd Dsuite
$ make
```

The Dsuite executable will be in the Build folder, so to run it type e.g. `./Build/Dsuite`; this will show the available commands. To execute e.g. the Dtrios command, type `./Build/Dsuite Dtrios`.

## [Optional] Installing the python3 Fbranch plotting script

If you want to plot the results of the f-branch calcuation (see below), you will need to install the python script for this using setuptools. You need an internet connection as some python dependencies will be downloaded from `pypi.org`. It may be necessary to exit python or conda virtual environments for this to work correctly.

```console
$ cd utils
$ python3 setup.py install --user --prefix=
```

The above should work on both mac and linux. Note that there is no text (not even whitespace) after the `=` above. If you want to use your own virtual environments, you can alternatively not run setup.py and just install the dependencies with `pip` or `conda`.


## Input files:
### Required files:
1. A [VCF](http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/) file, which can be compressed with gzip or bgzip. It can contain multiallelic loci and indels, but only biallelic loci will be used.
2. Population/species map `SETS.txt`: a text file with one individual per row and a tab separating the individual’s name from the name of the species/population it belongs to, as shown below:
```
Ind1    Species1
Ind2    Species1
Ind3    Species2
Ind4    Species2
Ind5    Species3
Ind6    Outgroup
Ind7    Outgroup
Ind8    xxx
...     ...
IndN    Species_n
```
At least one individual needs to be specified to be the outgroup by using the `Outgroup` keyword. If you want some individuals to be ignored entirely, use the `xxx` keyword instead. Therefore, you don't have to subset your VCF file if you only want to use a subset of the samples in it. 
### Optional files:
3. A tree in Newick format. The tree should have leaf labels corresponding to the species/population names. Branch lengths can be present but are not used.  
Valid examples:  
`(Species2,(Species1,(Species3,Species4)));`  
`(Species2:6.0,(Species1:5.0,(Species3:3.0,Species4:4.0)));`
4. The `test_trios.txt` file for `Dinvestigate`. One trio of populations/species per line, separated by a tab in the order `P1  P2  P3`:
```
Species1    Species2    Species3
Species1    Species4    Species2
...         ...         ...
```

## Commands (v0.3):
### Dsuite Dtrios - Calculate the D (ABBA-BABA) and f4-ratio (f_G) statistics for all possible trios of populations/species
```
Calculate the D (ABBA/BABA) and f4-ratio (f_G) statistics for all trios of species in the dataset (the outgroup being fixed)
the results are as definded in Patterson et al. 2012 (equivalent to Durand et al. 2011 when the Outgroup is fixed for the ancestral allele)
The SETS.txt should have two columns: SAMPLE_ID    SPECIES_ID
The outgroup (can be multiple samples) should be specified by using the keywork Outgroup in place of the SPECIES_ID

       -h, --help                              display this help and exit
       -k, --JKnum                             (default=20) the number of Jackknife blocks to divide the dataset into; should be at least 20 for the whole dataset
       -j, --JKwindow                          (default=NA) Jackknife block size in number of informative SNPs (as used in v0.2)
                                               when specified, this is used in place of the --JKnum option
       -r, --region=start,length               (optional) only process a subset of the VCF file
       -t, --tree=TREE_FILE.nwk                (optional) a file with a tree in the newick format specifying the relationships between populations/species
                                               D and f4-ratio values for trios arranged according to the tree will be output in a file with _tree.txt suffix
       -n, --run-name                          run-name will be included in the output file name

```
#### Output:
The output files with suffixes  `BBAA.txt`, `Dmin.txt`, and optionally `tree.txt` (if the `-t` option was used) contain the results: the D statistics, unadjusted p-values, and the f4-ratios (f_G). Please read the [manuscript](https://www.biorxiv.org/content/biorxiv/early/2019/05/10/634477.full.pdf) for more details. 

The output files with suffixes  `combine.txt` and  `combine_stderr.txt` are used as input to DtriosCombine. If you don't need to use DtriosCombine, you can safely delete these files.

### DtriosCombine - Combine results from Dtrios runs across genomic regions (e.g. per chromosome)
```
Combine the BBAA, ABBA, and BABA counts from multiple files (e.g per-chromosome) and output the overall D stats,
p-values and f4-ratio values

       -h, --help                              display this help and exit
       -n, --run-name                          run-name will be included in the output file name
       -t , --tree=TREE_FILE.nwk               (optional) a file with a tree in the newick format specifying the relationships between populations/species
                                               D and f4-ratio values for trios arranged according to the tree will be output in a file with _tree.txt suffix
       -s , --subset=start,length              (optional) only process a subset of the trios
```
#### Output:
As for `Dtrios`, there are output files with suffixes  `BBAA.txt`, `Dmin.txt`, and optionally `tree.txt` (if the `-t` option was used). They contain the overall combined results: the D statistics, unadjusted p-values, and the f4-ratios (f_G). 

###  Dinvestigate - Follow up analyses for trios with significantly elevated D: calculates D, f_d and f_dM in windows along the genome
```
Outputs D, f_d (Martin et al. 2014 MBE), and f_dM (Malinsky et al., 2015) in genomic windows
The SETS.txt file should have two columns: SAMPLE_ID    POPULATION_ID
The test_trios.txt should contain names of three populations for which the statistics will be calculated:
POP1   POP2    POP3
There can be multiple lines and then the program generates multiple ouput files, named like POP1_POP2_POP3_localFstats_SIZE_STEP.txt

       -h, --help                              display this help and exit
       -w SIZE,STEP --window=SIZE,STEP         (required) D, f_D, and f_dM statistics for windows containing SIZE useable SNPs, moving by STEP (default: 50,25)
       -n, --run-name                          run-name will be included in the output file name
```

###  Fbranch - A heuristic approach designed to aid the interpretation of many correlated f4-ratio results 
```
Usage: Dsuite Fbranch [OPTIONS] TREE_FILE.nwk FVALS_tree.txt
Implements the 'f-branch' type calculations developed by Hannes Svardal for Malinsky et al., 2018, Nat. Ecol. Evo.
Uses the f4-ratio (f_G) values produced by Dsuite Dtrios (or DtriosCombine) with the --tree option; this is the output of Dtrios with the "_tree.txt" suffix

       -h, --help                              display this help and exit
```
#### Output:
The f-branch statistic in matrix-like format. Use the plotting function below to display the f-branch statistic.

###  Plotting Fbranch
The output of `Dsuite Fbranch` can be plotted with `./utils/dtools.py` (see installation instructions above).

```
usage: dtools.py [-h] [-n RUN_NAME] [--outgroup OUTGROUP] [--use_distances]
                 [--ladderize]
                 fbranch tree

Plot f-branch statistic as produced by Dsuite. Produces .png and .svg files.

positional arguments:
  fbranch               Path to file containing f-branch matrix as produced by
                        Dsuite Fbranch.
  tree                  Path to .newick tree file as given to Dsuite Fbranch.

optional arguments:
  -h, --help            show this help message and exit
  -n RUN_NAME, --run-name RUN_NAME
                        Base file name for output plots. (default: fbranch)
  --outgroup OUTGROUP   Outgroup name in newick file. (default: Outgroup)
  --use_distances       Use actual node distances from newick file when
                        plotting tree. (default: False)
  --ladderize           Ladderize the input tree before plotting. (default:
                        False)
```

Running `dtools.py` yields a .png and an .svg file of the f-branch statistic along the input tree. You can edit the .svg file in a vector graphics editor (e.g., [inkscape](https://inkscape.org/) to your liking). See [Malinsky et al. 2018](https://www.nature.com/articles/s41559-018-0717-x) Fig. 3 and the Dsuite preprint biorxiv for examples and interpretation of f-branch plots.


## Change log:

```
Selected updates (full update history is accessible on gitHub):
v0.3 r21:   Automatic estimation of Jackknife window size to get a desired number of blocks
            Progress update in %
            f4-ratios are calculated by default by Dtrios 
            Updated documentation
v0.2 r20:   Subset option returns to DtriosCombine
v0.2 r19:   Fixed a bug in Fbranch where P1 and P2 positions where A in P2 positions and B in P1 positions were not considered
v0.2 r18:   Full implementation of D and f4-ratio in line with the Patterson et al. 2012 definitions 
            (affects only analyses where the outgroup allele is not fixed)
v0.2 r15:   Ironed bugs in Dinvestigate, truly useable from this point  
v0.2 r6:    First Fbranch version
v0.1 r1:    First workable Dsuite release 8th May 2019    

```
