# Population Genomics with R
----------------------
#### Author: Pierre Lesturgie
----------------------

## Analyse and filter genomics data using R functions
 
This repository provides python scripts to analyse and filter genomic data.
 
Installing the following packages will ensure all functions to work:
- parallel
- vcfR
- ggplot2
- ggthemes
- scales
- ade4
- rlist

## (1) Filtering and extracting info from VCFs
### Missing data
Returns a list with missing data per sample or SNPs or both
Need to input either a vcf or at GT matrix (obtained from vcfR extract.gt() function)

    missing.data(GT=NULL,vcf=NULL, SAMPLE=T, SNPs=T)

### Depth
Returns a list with depth of coverage data per sample or SNPs or both
Need to input either a vcf or a depth matrix (obtained from vcfR extract.gt() function)

    depth(DP=NULL,vcf=NULL, SAMPLE=T, SNPs=T)
    
### gc content
Returns GC data
Need to input either a vcf or at GT matrix (obtained from vcfR extract.gt() function)

    gc.content(GT=NULL, vcf=NULL)

### Heterozygosity
Returns a list with heterozygosity data per sample or SNPs or both
Need to input either a vcf or at GT matrix (obtained from vcfR extract.gt() function)

    heterozygosity(GT=NULL,vcf=NULL, SAMPLE=F, SNPs=T)

### Frequency of fixed alleles
Returns a list with frequency of fixed alleles data per sample or SNPs or both
Need to input either a vcf or at GT matrix (obtained from vcfR extract.gt() function)

    fixed.alleles(GT=NULL,vcf=NULL, SAMPLE=F, SNPs=T)

### Minor allele frequency
Returns MAF data
Need to input either a vcf or at GT matrix (obtained from vcfR extract.gt() function)

    minor.allele.frequency(GT=NULL,vcf=NULL)

### Pick random SNP per locus
Returns a filtered VCF. This is usually meaningful for RAD-like data (i.e., where there are many loci)

    random.snp(vcf)

## (2) Site frequency spectrum (SFS)

### Folded SFS calculation
Returns folded SFS (accross all SNPs or in sliding windows)
Need to input either a vcf or at GT matrix (obtained from vcfR extract.gt() function)
If needing monomorphinc sites (i.e., with a minor allele frequency of 0) use <mono=T>
Possibility of parallelization 

    sfs.folded(vcf=NULL,GT=NULL,slide=100000,jump=25000,sliding_window=FALSE,paral=F,mono=F)

### Unfolded SFS calculation
Returns folded SFS (and possibly a derived vcf)
Need to input either a vcf or at GT matrix (obtained from vcfR extract.gt() function) and the sample name of the outgroup

    sfs.unfolded(vcf=NULL, GT=NULL, outgroup, return_derived_vcf=F)

### SFS folding 
Returns folded SFS from unfolded

    fold.SFS(sfs)

### Calculates a correct SFS with missing data
Returns a folded SFS from a vcf and a maximum missing data rate

    corr.sfs.missing.data(vcf,md)

### SFS normalization as in Lapierre et al. (2021)
Returns normalized SFS from sfs (Without monomorphic sites)
Possibility to return the expected norm SFS and expected SFS under and panmictic and constant population scenario

    normalized.expected.SFS(sfs, return.expected.norm=TRUE, return.norm=TRUE, return.expected.SFS=TRUE)

### SFS normalization as in Lapierre et al. (2021) in sliding windows
Returns normalized SFS from sfs (Without monomorphic sites) per windows

    normalized.expected.SFS.sliding.windows(sfs)

### Euclidian distance in sliding windows
Returns Euclidian distance between X and Y per windows (typically two SFS)

    euclidian.distance.sliding.window(X,Y=avgSFS)
    
### Distance SFS in windows to Average SFS 
Input is a data frame with the first three columns being: <CONTIG>, <low_bound_window>, <high_bound_window>
followed by the SFS for each window

    dist.normSFS(data)

### Calculates Watterson's (1979) estimate of genetic diversity
Input is simply a folded SFS (only polymorphic sites)

    watterson(sfs)


### Psi calculation as in Peter & Slatkin (2013, 2015)
Input is a dataframe with number of derived allees per SNP and individual. 
The second script (psi2dist) returns a distance matrix with PSI valyes and significancy

    psi(der_all, bootstrap=1000)
    psi2dist(psi,signif=T)


### Genetic distance between two individuals
There are two functions: 

#### (A) : D = 1 - (shared_alleles / 2)
Input is simply a genotype table (obtained using vcfR extract.gt function)

    genetic.distance(gt)

#### (B) : Bray-Curtis distance
Input is simply a vcf

    BC.distance.individual(vcf,bootstrap=NULL)

### Plot Stairway Plot with confidence intervals
Input is a **named** list of multiple output dataframe from stairwayplot (*.summary)
Also works with a single dataframe but still needs to be stored in a list. 

    plot.stairway.IC(data_list,var="Ne",output='stairway.pdf',alpha=0.1,cols=NULL,CI=T,xlim=NULL, legend=T,
                           leg_pos=c(0.78, 0.92),ylim=NULL,x.breaks=4,y.breaks=5,ncol.leg=1,by.gen=NULL,
                           vline=NULL,col_vline='grey',ltyp=NULL,ltyp_vline='dashed')

### Neighbouring SFS in ABC sumstat
Computes the N closest SFS to the observed one in one or multiple simulated models 

    neighbours.SFS(sumstat,target,n.closest=100,plot.output=NULL)
    multi.neigh.SFS(sumstat_list,target,n.closest=100,plot.output=NULL)

### Fluctuations in stairway plot output
computes the sum of slopes between discretized intervals from the stairway plot output

    fluctuations.stairway(stairway,N_DISCRETIZE)


## (3) FST
### Computes Hudson's (19XX) pairwise-FST
This returns:
- Overall pairwise FST between all locations
- N resampling values for each pairwise comparison (optional)
- Nucleotide pairwise FST between all locations
- Pairwise FST values in sliding windows (optional)

Input is a VCF and a list of populations

    fst.hudson(vcf, pop_list,resampling=100,sliding_window=FALSE,slide=NULL,jump=NULL,write=FALSE)

### Computes Isolation By Distance using a Mantel Test
Input are a geographic distance matrix and a FST distance matrix

    ibd.fst(dgeo,dfst)

##### The first one can be obtained from geo coordinates (lon, lat) with:

    dgeo = distance.geo(map, unit="km")
    
##### The second inputing pairwise FST in: 

    dfst = dist.fst(df_fst, popmap=NULL, signif=F, pv=0.01, sep_pop="/")

### Plot heatmap and similarity graphy of FST
Input is a distance matrix of FST

    plot.fst(FST,output="",gradmid=F)


