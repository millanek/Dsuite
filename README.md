#  Dsuite
Publication:  
Malinsky, M., Matschiner, M. and Svardal, H. (2021) Dsuite ‐ fast D‐statistics and related admixture evidence from VCF files. Molecular Ecology Resources 21, 584–595. doi: [https://doi.org/10.1111/1755-0998.13265](https://doi.org/10.1111/1755-0998.13265)  
Free to view author [link](https://onlinelibrary.wiley.com/share/author/QNEE6JI7DGUSBA4Y8ZGU?target=10.1111/1755-0998.13265)  
Malawi cichlid data used in the manuscript:
[VCF file](https://github.com/millanek/tutorials/blob/master/DsuiteData/Malinsky_et_al_2018_LakeMalawiCichlids_scaffold_0.vcf.gz); [SETS.txt file](https://github.com/millanek/tutorials/blob/master/DsuiteData/sets.txt)<br>
Simulated 20-species data used in the manuscript: [VCF file](https://github.com/millanek/tutorials/blob/master/DsuiteData/with_geneflow.vcf.gz); [SETS.txt file](https://github.com/millanek/tutorials/blob/master/DsuiteData/species_sets.txt); [TREE_FILE.nwk](https://github.com/millanek/tutorials/blob/master/DsuiteData/simulated_tree_with_geneflow.nwk) (input tree)

There is also a very detailed [tutorial](https://github.com/millanek/tutorials/tree/master/analysis_of_introgression_with_snp_data) that I prepared with input from [@mmatschiner](https://github.com/mmatschiner).

## Quickstart:
```
Commands:
           Dtrios                  Calculate D (ABBA-BABA) and f4-ratio statistics for all possible trios of populations/species
           DtriosCombine           Combine results from Dtrios runs across genomic regions (e.g. per-chromosome)
           Dinvestigate            Follow up analyses for trios with significantly elevated D:
                                   calculates f_d, f_dM, and d_f in windows along the genome
           Fbranch                 Calculate D and f statistics for branches on a tree that relates the populations/species

Experimental:
           Dquartets               Calculate D (ABBA-BABA) and f4-ratio statistics for all possible quartets of populations/species
                                   (no outgroup specified)

Usage:
a) Dsuite Dtrios [OPTIONS] INPUT_FILE.vcf SETS.txt
b) Dsuite Dquartets [OPTIONS] INPUT_FILE.vcf SETS.txt
c) Dsuite Dinvestigate [OPTIONS] INPUT_FILE.vcf.gz SETS.txt test_trios.txt
d) Dsuite Fbranch [OPTIONS] TREE_FILE.nwk FVALS_tree.txt
```

## Input files:
### Required files:
1. A [VCF](http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/) file, which can be compressed with gzip or bgzip. It can contain multiallelic loci and indels, but only biallelic SNPs will be used.
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
If you want some individuals to be ignored, use the `xxx` keyword. Therefore, you don't have to subset your VCF file if you want to use only a subset of the samples in it.

For `Dtrios`, at least one individual needs to be specified to be the outgroup by using the `Outgroup` keyword as shown above.

All species/populations are treated equally in `Dquartets` - there should not be any outgroup.
### Optional files:
3. A tree in Newick format. The tree should have leaf labels corresponding to the species/population names. Branch lengths can be present but are not used. To use  `Fbranch`, the tree must be rooted using the Outgroup. 
Valid examples:  
`(Species2,(Species1,(Species3,Species4)));`  
`(Species2:6.0,(Species1:5.0,(Species3:3.0,Species4:4.0)));`
4. The `test_trios.txt` file for `Dinvestigate`. One trio of populations/species per line, separated by a tab in the order `P1  P2  P3`:
```
Species1    Species2    Species3
Species1    Species4    Species2
...         ...         ...
```
### Piped VCF input:
It is possible to 'pipe' the genotype data into  `Dsuite Dtrios` or `Dsuite Dquartets` from another program, such as bcftools, allowing custom filtering of the VCF file.  Just use the `stdin` keyword in place of the VCF file name. It is also necessary to provide the number of lines in the filtered VCF via the  `-l` option to the Dsuite programs. For example, to filter a VCF for overall mimimum depth of at least 1000 across all samples, you would use the following commands:
```
NUMLINES=$(bcftools view -i 'INFO/DP>1000' INPUT_FILE.vcf.gz | wc -l)  # to get NUMLINES
bcftools view -i 'INFO/DP>1000' INPUT_FILE.vcf.gz | Dsuite Dtrios -l $NUMLINES stdin SETS.txt
```

## Installation
### Main program:
To compile you must have a reasonably recent GCC (>=4.9.0) or clang compiler (on mac OS this comes with Command Line Tools) and the zlib compression library (https://www.zlib.net). Both should already be present on most systems. 

```console
$ git clone https://github.com/millanek/Dsuite.git
$ cd Dsuite
$ make
```

The Dsuite executable will be in the Build folder, so to run it type e.g. `./Build/Dsuite`; this will show the available commands. To execute e.g. the Dtrios command, type `./Build/Dsuite Dtrios`.

### [Optional] Installing the python3 Fbranch plotting script

If you want to plot the results of the f-branch calcuation (see below), you will need to install the python script for this using setuptools. You need an internet connection as some python dependencies will be downloaded from `pypi.org`. It may be necessary to exit python or conda virtual environments for this to work correctly.

```console
$ cd utils
$ python3 setup.py install --user --prefix=
```

The above should work on both mac and linux. Note that there is no text (not even whitespace) after the `=` above. If you want to use your own virtual environments, you can alternatively not run setup.py and just install the dependencies with `pip` or `conda`.


## Commands (v0.5 r44):
### Dsuite Dtrios - Calculate the D (ABBA-BABA) and f4-ratio statistics for all possible trios of populations/species
```
Usage: Dsuite Dtrios [OPTIONS] INPUT_FILE.vcf SETS.txt
Calculate the D (ABBA/BABA) and f4-ratio statistics for all trios of species in the dataset (the outgroup being fixed)
The results are as definded in Patterson et al. 2012 (equivalent to Durand et al. 2011 when the Outgroup is fixed for the ancestral allele)
The SETS.txt should have two columns: SAMPLE_ID    SPECIES_ID
The outgroup (can be multiple samples) should be specified by using the keywork Outgroup in place of the SPECIES_ID


      -h, --help                              display this help and exit
      -k, --JKnum                             (default=20) the number of Jackknife blocks to divide the dataset into; should be at least 20 for the whole dataset
      -j, --JKwindow                          (default=NA) Jackknife block size in number of informative SNPs (as used in v0.2)
                                              when specified, this is used in place of the --JKnum option
      -r, --region=start,length               (optional) only process a subset of the VCF file; both "start" and "length" indicate variant numbers
                                              e.g. --region=20001,10000 will process variants from 20001 to 30000
      -t, --tree=TREE_FILE.nwk                (optional) a file with a tree in the newick format specifying the relationships between populations/species
                                              D and f4-ratio values for trios arranged according to the tree will be output in a file with _tree.txt suffix
      -o, --out-prefix=OUT_FILE_PREFIX        (optional) the prefix for the files where the results should be written
                                              output will be put in OUT_FILE_PREFIX_BBAA.txt, OUT_FILE_PREFIX_Dmin.txt, OUT_FILE_PREFIX_tree.txt etc.
                                              by default, the prefix is taken from the name of the SETS.txt file
      -n, --run-name                          (optional) run-name will be included in the output file name after the PREFIX
      --no-f4-ratio                           (optional) don't calculate the f4-ratio
      -l NUMLINES                             (optional) the number of lines in the VCF input - required if reading the VCF via a unix pipe
      -g, --use-genotype-probabilities        (optional) use probabilities (GP tag) or calculate them from likelihoods (GL or PL tags) using a Hardy-Weinberg prior
                                              the probabilities are used to estimate allele frequencies in each population/species
      -p, --pool-seq=MIN_DEPTH                (optional) VCF contains pool-seq data; i.e., each 'individual' is a population
                                              allele frequencies are then estimated from the AD (Allelic Depth) field, as long as there are MIN_DEPTH reads
                                              e.g MIN_DEPTH=5 may be reasonable; when there are fewer reads, the allele frequency is set to missing
      -c, --no-combine                        (optional) do not output the "_combine.txt" and "_combine_stderr.txt" files
                                              these are needed only for DtriosCombine
      --KS-test-for-homoplasy                 (optional) Test whether strong ABBA-informative sites cluster along the genome
```
#### Output:
The output files with suffixes  `BBAA.txt`, `Dmin.txt`, and optionally `tree.txt` (if the `-t` option was used) contain the results: the D statistics, Zscore, unadjusted p-values, the f4-ratios, and counts of the BBAA, BABA, and ABBA patterns. Please read the [manuscript](https://doi.org/10.1111/1755-0998.13265) for more details. 

The output files with suffixes  `combine.txt` and  `combine_stderr.txt` are used as input to DtriosCombine. If you don't need to use DtriosCombine, you can safely delete these files.

### DtriosParallel

We provide a python script for parallel execution at `<Dsuite_path>/utils/DtriosParallel`. The usage is analogous to Dsuite Dtrios, except that the order of `SETS.txt` and `INPUT_FILE.vcf` is swapped in the command line so that the user can optionally provide multiple VCF files (whitespace separated). The script will autmatically call `DtriosCombine` to combine results of all VCF files into a single set of results files.

```
$ ./utils/DtriosParallel --help
usage: DtriosParallel [-h] [-k JKNUM] [-j JKWINDOW] [-t TREE] [-n RUN_NAME]
                      [-l NUMLINES] [-g] [-p --pool-seq=MIN_DEPTH] [-c]
                      [--cores CORES] [--keep-intermediate]
                      [--logging_level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                      [--dsuite-path DSUITE_PATH]
                      [--environment-setup ENVIRONMENT_SETUP]
                      SETS.txt INPUT_FILE.vcf [INPUT_FILE.vcf ...]

This python script automates parallelisation of Dsuite Dtrios/ Dsuite
DtriosCombine. The usage is analogous to Dsuite Dtrios but computation is
performed on multiple cores (default: number of available CPUs). ATTENTION:
The order of SETS.txt and INPUT_FILE.vcf is swapped compared to Dsuite Dtrios.
This is so that multiple VCF input files can be provided. All vcf files should
have the same samples, (e.g., different chromosomes of the same callset).
Output_files are placed in the same folder as as the the SETS.txt file and
named DTParallel_<SETS_basename>_<run_name>_combined_BBAA.txt etc. This script
should run on most systems with a standard python installation (tested with
python 2.7 and 3.6).


positional arguments:
  SETS.txt              The SETS.txt should have two columns: SAMPLE_ID
                        SPECIES_ID The outgroup (can be multiple samples)
                        should be specified by using the keyword Outgroup in
                        place of the SPECIES_ID
  INPUT_FILE.vcf        One or more whitespace separated SNP vcf files.

optional arguments:
  -h, --help            show this help message and exit
  -k JKNUM, --JKnum JKNUM
                        (default=20) the number of Jackknife blocks to divide
                        the dataset into; should be at least 20 for the whole
                        dataset
  -j JKWINDOW, --JKwindow JKWINDOW
                        Jackknife block size in number of informative SNPs (as
                        used in v0.2) when specified, this is used in place of
                        the --JKnum option
  -t TREE, --tree TREE  a file with a tree in the newick format specifying the
                        relationships between populations/species D and
                        f4-ratio values for trios arranged according to the
                        tree will be output in a file with _tree.txt suffix
  -n RUN_NAME, --run-name RUN_NAME
                        run-name will be included in the output file name
  -l NUMLINES           (optional) the number of lines (SNPs) in the VCF
                        input(s) - speeds up operation if known. If N
                        INPUT_FILE.vcf files are provided, there must be N
                        comma-separated integers provided without whitespace
                        between them.
  -g, --use-genotype-probabilities
                        (optional) use probabilities (GP tag) or calculate
                        them from likelihoods (GL or PL tags) using a Hardy-
                        Weinberg prior
  -p --pool-seq=MIN_DEPTH, --pool-seq --pool-seq=MIN_DEPTH
                        (default=20) the number of Jackknife blocks to divide
                        the dataset into; should be at least 20 for the whole
                        dataset
  -c, --no-combine      (optional) do not run DtriosCombine to obtain a single
                        combined results file
  --cores CORES         (default=CPU count) Number of Dsuite Dtrios processes
                        run in parallel.
  --keep-intermediate   Keep region-wise Dsuite Dtrios results.
  --logging_level {DEBUG,INFO,WARNING,ERROR,CRITICAL}, -v {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Minimun level of logging.
  --dsuite-path DSUITE_PATH
                        Explicitly set the path to the directory in which
                        Dsuite is located. By default the script will first
                        check whether Dsuite is accessible from $PATH. If not
                        it will try to locate Dsuite at ../Build/Dsuite.
  --environment-setup ENVIRONMENT_SETUP
                        Command that should be run to setup the environment
                        for Dsuite. E.g., 'module load GCC' or 'conda


```

#### Output:

Output are \_BBAA, \_Dmin files analogus to Dtrios/DtriosCombine and are placed in the same folder as as the the `SETS.txt` file and
named `DTParallel_<SETS_basename>_<run_name>_combined_BBAA.txt` etc.

### DtriosCombine - Combine results from Dtrios runs across genomic regions (e.g. per chromosome)
```
Usage: Dsuite DtriosCombine [OPTIONS] DminFile1 DminFile2 DminFile3 ....

Combine the BBAA, ABBA, and BABA counts from multiple files (e.g per-chromosome) and output the overall D stats,
p-values and f4-ratio values

       -h, --help                              display this help and exit
       -o, --out-prefix=OUT_FILE_PREFIX        (optional) the prefix for the files where the results should be written
                                               output will be put in OUT_FILE_PREFIX_combined_BBAA.txt, OUT_FILE_PREFIX_combined_Dmin.txt, OUT_FILE_PREFIX_combined_tree.txt etc.
                                               by default, the prefix is "out"
       -n, --run-name                          (optional) run-name will be included in the output file name after the PREFIX
       -t , --tree=TREE_FILE.nwk               (optional) a file with a tree in the newick format specifying the relationships between populations/species
                                               D and f4-ratio values for trios arranged according to the tree will be output in a file with _tree.txt suffix
       -s , --subset=start,length              (optional) only process a subset of the trios
```
#### Output:
As for `Dtrios`, there are output files with suffixes  `BBAA.txt`, `Dmin.txt`, and optionally `tree.txt` (if the `-t` option was used). They contain the overall combined results: the D statistics, Zscore, unadjusted p-values, the f4-ratios, and counts of the BBAA, BABA, and ABBA patterns.

###  Dinvestigate - Follow up analyses for trios with significantly elevated D: calculates D, f_d and f_dM in windows along the genome
```
Usage: Dsuite Dinvestigate [OPTIONS] INPUT_FILE.vcf.gz SETS.txt test_trios.txt

Outputs D, f_d (Martin et al. 2014 MBE), f_dM (Malinsky et al., 2015), and d_f (Pfeifer & Kapan, 2019) in genomic windows
The SETS.txt file should have two columns: SAMPLE_ID    POPULATION_ID
The test_trios.txt should contain names of three populations for which the statistics will be calculated:
POP1   POP2    POP3
There can be multiple lines and then the program generates multiple ouput files, named like POP1_POP2_POP3_localFstats_SIZE_STEP.txt

       -h, --help                              display this help and exit
       -w SIZE,STEP --window=SIZE,STEP         (required) D, f_D, f_dM, and d_f statistics for windows containing SIZE useable SNPs, moving by STEP (default: 50,25)
       -n, --run-name                          run-name will be included in the output file name
```

###  Fbranch - A heuristic approach designed to aid the interpretation of many correlated f4-ratio results 
```
Usage: Dsuite Fbranch [OPTIONS] TREE_FILE.nwk FVALS_tree.txt
Implements the 'f-branch' type calculations developed by Hannes Svardal for Malinsky et al., 2018, Nat. Ecol. Evo.
Uses the f4-ratio (f_G) values produced by Dsuite Dtrios (or DtriosCombine) with the --tree option; this is the output of Dtrios with the "_tree.txt" suffix
To use  Fbranch, the tree in TREE_FILE.nwk must be rooted with the Outgroup.
Output to stdout

      -p, --pthresh                           (default=0.01) fb scores whose associated p-value is less than 
      -Z, --Zb-matrix                         (optional)  output the equivalent of fb-statistic, but with Z-scores to assess statistical significance
                                              this will be printed below the f-branch matrix
      -h, --help                              display this help and exit
```

#### Output:
The f-branch statistic in matrix-like format. Use the plotting function below to display the f-branch statistic.

###  Plotting Fbranch
The output of `Dsuite Fbranch` can be plotted with `./utils/dtools.py` (see installation instructions above).

```
usage: dtools.py [-h] [-n RUN_NAME] [--outgroup OUTGROUP] [--use_distances]
                 [--ladderize]
                 fbranch.txt tree.newick

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

Running `dtools.py` yields a .png and an .svg file of the f-branch statistic along the input tree. You can edit the .svg file in a vector graphics editor (e.g., [inkscape](https://inkscape.org/) to your liking). See [Malinsky et al. 2018](https://www.nature.com/articles/s41559-018-0717-x) Fig. 3 and the Dsuite [paper](https://doi.org/10.1111/1755-0998.13265) for examples and interpretation of f-branch plots.

### (experimental) Dsuite Dquartets - Calculate the D (ABBA-BABA) and f4-ratio statistics for all possible quartets of populations/species (no outgroup)
```
Usage: Dsuite Dquartets [OPTIONS] INPUT_FILE.vcf SETS.txt
Calculate the D (ABBA/BABA) and f4-ratio (f_G) statistics for all quartets of species in the dataset (there is no outgroup)
The results are as definded in Patterson et al. 2012
The SETS.txt should have two columns: SAMPLE_ID    SPECIES_ID

-h, --help                              display this help and exit
-k, --JKnum                             (default=20) the number of Jackknife blocks to divide the dataset into; should be at least 20 for the whole dataset
-j, --JKwindow                          (default=NA) Jackknife block size in number of informative SNPs (as used in v0.2)
                                        when specified, this is used in place of the --JKnum option
-r, --region=start,length               (optional) only process a subset of the VCF file
-t, --tree=TREE_FILE.nwk                (optional) a file with a tree in the newick format specifying the relationships between populations/species
                                        D and f4-ratio values for trios arranged according to the tree will be output in a file with _tree.txt suffix
-o, --out-prefix=OUT_FILE_PREFIX        (optional) the prefix for the files where the results should be written
                                        output will be put in OUT_FILE_PREFIX_BBAA.txt, OUT_FILE_PREFIX_Dmin.txt, OUT_FILE_PREFIX_tree.txt etc.
                                        by default, the prefix is taken from the name of the SETS.txt file
-n, --run-name                          (optional; default=quartets) run-name will be included in the output file name after the PREFIX
--no-f4-ratio                           (optional) don't calculate the f4-ratio
-l NUMLINES                             (optional) the number of lines in the VCF input - required if reading the VCF via a unix pipe
```

### Parallelisation with DtriosParallel

This python script, included in the `./utils/` subfolder, automates parallel runs of `Dsuite Dtrios` across multiple cores on one computer and automatically combines the results using `Dsuite DtriosCombine`. The usage is analogous to `Dsuite Dtrios` (currently with more limited options) but computation is performed on multiple cores (default: number of available CPUs). It should run on most systems with a standard python installation (tested with python 2.7 and 3.6).

```
DtriosParallel [-h] [--cores CORES] [-k JKNUM] [-j JKWINDOW] [-t TREE]
                      [-n RUN_NAME] [--keep-intermediate]
                      [--logging_level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                      [--dsuite-path DSUITE_PATH]
                      [--environment-setup ENVIRONMENT_SETUP]
                      INPUT_FILE.vcf SETS.txt


positional arguments:
  INPUT_FILE.vcf
  SETS.txt              The SETS.txt should have two columns: SAMPLE_ID
                        SPECIES_ID The outgroup (can be multiple samples)
                        should be specified by using the keyword Outgroup in
                        place of the SPECIES_ID

optional arguments:
  -h, --help            show this help message and exit
  --cores CORES         (default=CPU count) Number of Dsuite Dtrios processes
                        run in parallel.
  -k JKNUM, --JKnum JKNUM
                        (default=20) the number of Jackknife blocks to divide
                        the dataset into; should be at least 20 for the whole
                        dataset
  -j JKWINDOW, --JKwindow JKWINDOW
                        Jackknife block size in number of informative SNPs (as
                        used in v0.2) when specified, this is used in place of
                        the --JKnum option
  -t TREE, --tree TREE  a file with a tree in the newick format specifying the
                        relationships between populations/species D and
                        f4-ratio values for trios arranged according to the
                        tree will be output in a file with _tree.txt suffix
  -n RUN_NAME, --run-name RUN_NAME
                        run-name will be included in the output file name
  --keep-intermediate   Keep region-wise Dsuite Dtrios results.
  --logging_level {DEBUG,INFO,WARNING,ERROR,CRITICAL}, -l {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Minimun level of logging.
  --dsuite-path DSUITE_PATH
                        Explicitly set the path to the directory in which
                        Dsuite is located. By default the script will first
                        check whether Dsuite is accessible from $PATH. If not
                        it will try to locate Dsuite at ../Build/Dsuite.
  --environment-setup ENVIRONMENT_SETUP
                        Command that should be run to setup the environment
                        for Dsuite. E.g., 'module load GCC' or 'conda
                        activate'
```

## Change log:

```
Selected updates (full update history is accessible on gitHub):
v0.5 r48:   --KS-test-for-homoplasy now works accurately for sample sizes >= 16,000  
v0.5 r47:   --KS-test-for-homoplasy output now has more accurate p-values (one-sample KS-test for uniformity)
v0.5 r46:   Support for arbitrary ploidy in Dtrios
v0.5 r45:   BUG FIX for r44 where the trio orientation in the "_tree.txt" output files was wrong (P1 and P2 swapped). This is fixed now.
v0.5 r44:   Major update:   - code re-factoring, including proper subsampling for f4-ratio calculations
                            - First implementation of the Kolgomorov-Sminov test for homoplasy (--KS-test-for-homoplasy) in Dtrios; still somewhat experimental and works only in the "_BBAA.txt" output
                            - Very small p-values now don't get rounded to 0 but are bounded at 2.3e-16
v0.4 r43:   First implementation of the pool-seq (-p) option in in Dtrios 
v0.4 r28:   Merged DtriosParallel from https://github.com/feilchenfeldt and refreshed documentation
v0.3 r27:   Added the -o (--out-prefix) option to allow more flexibility in naming output files
v0.3 r25:   Added the Dquartets program - D and f4-ratio calculation without any outgroup, for all quartets of populations/species
v0.3 r24:   Z-scores and site pattern counts (BBAA, ABBA, BABA) are now in the output of Dtrios, DtriosCombine, and Dquartets 
v0.3 r23:   Allow piped (stdin) VCF format input to Dtrios; this facilitates e.g. pre-filtering and/or bcf input using bcftools 
v0.3 r22:   Adding the d_f statistic (Pfeifer & Kapan, 2019) to Dinvestigate
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
