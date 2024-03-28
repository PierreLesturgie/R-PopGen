# Population Genomics with R
----------------------
#### Author: Pierre Lesturgie
----------------------

## Analyse and filter genomics data using R functions
 
 This repository provides python scripts to analyse and filter genomic data from VCF files. 
 
 For each script, all arguments are indicated using the --help argument.  
 
 Scripts run with python 3.8

### (1) Filtering SNPs with excess heterozygotes
Example with a threshold at 80% (output is  a vcf file named output.vcf)

	python HET.py --input input.vcf.gz -H 0.8 -O output 

### (2) Binning the vcf (for linkage disequilibrium)
There are two ways of binning: 
#### 1 - By regions, i.e., only keeps regions of <region> bp distant of <bin> from each other
Example with with regions of 100 bp distant of 100,000 bp to each other

	python binning.py --input input.vcf.gz --output binned_100k_regions_100bp.vcf.gz --bin 100000 --region 100

#### 2 - By SNP, i.e., only keeps SNPs distant of <bin> from each other
Example with with SNPs distant of 100,000 bp to each other

	python binning.py --input input.vcf.gz --output binned_100k.vcf.gz --bin 100000

### (3) Polarize alleles of a vcf for a given outgroup individual
Example for a polarization of alleles for outgroup named "GN18033" (outputs: polarized.vcf.gz and polarized.der)

	python Polarize.py --input input_with_outgroup.vcf.gz --outgroup GN18033 --output polarized -t tempola -D

### (4) Site Frequency Spectrum 
There are two ways of computing the SFS: 
#### 1 - Accross all sites 
Remove --folded for unfolded sfs

	python SFS.py --input input.vcf.gz --folded --output folded.sfs

#### 2 - In sliding windows
Remove --folded for unfolded sfs
Example for contig SUPER_1 (--chr SUPER_1) in windows of 1,000,000 bp with a jump of 500,000 bp

	python SFS.py --input input.vcf.gz --folded --output sliding_folded.sfs --chr SUPER_1 --sliding_window --window 1000000 --jump 500000
