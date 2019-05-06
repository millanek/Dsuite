#  Dsuite
Fast calculation of D-statistics and related admixture evidence directly from VCF files.


## REQUIREMENTS

Assuming that you are in a Linux/UNIX like environment, you need:
1. an up-to-date version of the gcc compiler, e.g. "gcc (Ubuntu 4.4.3-4ubuntu5) 4.4.3"
2. GNU Scientific Development Library, (libgsl0-dev and libgsl0 in Ubuntu, gsl-devel in OpenSUSE) 
3. Mac/Apple users need to have the Command Line Tools installed

## Compilation

To compile you must have a resonably recent GCC or clang compiler (on mac OS this comes with Command Line Tools) and the zlib compression library (https://www.zlib.net). Both should already be present on most systems. 


Then simply change directories into the Dsuite folder, and type:
`make`
The Dsuite executable will be in the Build folder.

## Input files:
1. Bialllelic variants in a VCF file. The fille can be gzipped.
http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/
2. 


