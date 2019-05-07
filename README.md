#  Dsuite
Fast calculation of D-statistics and related admixture evidence directly from VCF files.

## Compilation

To compile you must have a reasonably recent GCC or clang compiler (on mac OS this comes with Command Line Tools) and the zlib compression library (https://www.zlib.net). Both should already be present on most systems. 


Then simply change directories into the Dsuite folder, and type:
`make`

The Dsuite executable will be in the Build folder.

## Input files:
### Required files:
1. Bialllelic variants in a VCF file. The file can be gzipped.
http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/
2. Population/species map `SETS.txt`: a text file with one individual per row and a tab separating the individual’s name from the name of the species/population it belongs to, as shown below:
```
Ind1    Species1
Ind2    Species1
Ind3    Species2
Ind4    Species2
Ind5    Species3
...     ...
IndN    Species_n
```
### Optional files:
3. A tree in Newick format. For now I only have a very basic Newick format parser. Therefore, the tree should only have leaf labels and no branch lengths.
4. The `test_trios.txt` file for `Dinvestigate`. One trio of populations/species per line, separated by a tab in the order `P1  P2  P3`:
```
Species1    Species2    Species3
Species1    Species4    Species2
...         ...         ...
```



