# WISH

Weighted Interaction SNP Hub R Package

# Installation

Install WISH with the following commands:

```
install.packages(c("devtools","curl", "httr"))

library("devtools")

install_github("AQS-Group/WISH")
```
For the reference manual see the WISH.pdf file

# Quick Start

For using WISH you need to have your genotype data in the plink format.
See here for information:

http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped

Once you have a ped file you must make sure that is is 0,1,2 coded.
This can be done using the plink recode function:

plink --file \<input ped file\> --recode12

Further you need a transposed ped file. This is create with the following
command:

plink --file \<input 0,1,2 coded ped file\> --recode --transpose

The next step is to generat a genotype matrix in R using thegenerat.genotype() function. 
This will transform the alleles into single genotypes. Note that at this stage 
you should give a list of selected SNPs IDs or associated p-values for filtering
if you have not previously filtered your data (see the generate.genotype() documentation for
details on filtering. We would recomend to somewhat maximise the number of interactions calculated, 
with 10.000-20.000 SNPs being fairly easy with to run with ~10-15 cores, and more is possible
with strong computing facilities.


```
genotype <-generate.genotype(<input ped file path>,<input tped file path>,gwas.id=<selected list of id>,gwas.p=<p-values of input SNPs>)
```

***warning*** If you have more than about 1 million SNPs you must either fast.read = F which will slow down the loading time significantly.  You can also increase your stacklimit using ulimit in the command line,but do this only if you know what you are doing. 

After generating the genotypes file it is recomended to run a test run to estimate run time
of the epistatic interaction calculation based on available computing setup:
```
epistatic.correlation(<phenotype dataframe>, genotype,parallel = <number of cores available> ,test=T)
```
This will give you a indication of expected run time given your input. The next step is to run the analysis:
We recommend using simple=F for better results:
```
correlations<-epistatic.correlation(<phenotype dataframe>, genotype,parallel = <number of cores available> ,test=F,simple=F)
```
Once you have calculated epistatic correlations you can get a course grained overview of the results using
the genome.interaction() function:
```
genome.interaction(<input tped file>, correlations)
```
If you want to compare individual chromosomes see the pairwise.chr.map() function

Finally to create modules of interacting markers use the generate.modules function:
```
modules <-generate.modules(correlations)
```
***warning*** If you have epistatic correlatoin coefficients values that are strong outliers this can heavily affect
this step or even make it fail. Thus it is recommended to set those coefficients to 0, for example if we only want 
coefficients between 1000 and -1000:
```
correlations$Coefficients[correlations$Coefficients > 1000 |  correlations$Coefficients < -1000] <-0
```
The range of desired values is individual depending on the properties of each dataset.

# Contact
Victor A. O Carmelo, vaoc@sund.ku.dk
Faculty of Health and Medical Sciences
Department of Veterinary and Animal Sciences
Animal Breeding, Quantitative Genetics & Systems Biology Group

# References

Lisette J.A. Kogelman and Haja N.Kadarmideen (2014). 
Weighted Interaction SNP Hub (WISH) network method for building genetic 
Networks for complex diseases and traits using whole genome genotype data. 
BMC Systems Biology 8(Suppl 2):S5.  
http://www.biomedcentral.com/1752-0509/8/S2/S5.