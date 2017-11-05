# WISH

Weighted Interaction SNP Hub R Package

# Installation

Install WISH with the following commands:

For All users Users:

```
### Installing WGCNA first
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("GO.db", "preprocessCore", "impute"))
install.packages("WGCNA")
### Instaling devtools
install.packages(c("devtools","curl", "httr"))
### Install WISH
library("devtools")
install_github("QSG-Group/WISH")
```

If there are problems try this way:

```
### Instaling devtools
install.packages(c("devtools","curl", "httr"))


### Installing WGCNA first
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("GO.db", "preprocessCore", "impute"))
install.packages("WGCNA")

### Install Rest of Dependencies
install.packages(c("doParallel", "foreach","fastcluster", "Rcpp", "RcppEigen", "data.table", "corrplot", "heatmap3", "flashClust", "bigmemory", "parallel"))

### Install WISH
library("devtools")
install_github("QSG-Group/WISH")
```
For the reference manual see the WISH.pdf file
# Test files

The files test.ped, test.tped and test_pheno.txt show how you have to structure your input files
and allow you to test the comands.

# Quick Start
For using WISH you need to have your genotype data in the plink format.
See here for information on the plink data format:

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

***Data Pre-filtering***

We recommend prefiltering your data using a main effect filter. For example you can run a simple GWAS using plink:

plink --file \<ped file basename\> --linear --o \<output basename\>

The computed p-values can be used in R as filter. 

***Loading data into R***

The functions in the WISH R package accept both filepaths or data frames as input. To load the data easily into R
use following commands:
```
library(data.table)
ped <- tped <- fread(<filepath to pedfile>, data.table = F)
tped <- fread(<filepath to tpedfile>, data.table = F)
```
If memory load is a problem it is recommended to use the file paths, as the working enviroment
is duplicated when using multiple threads.
There is no strict guideline for the p-value threshold, but 
given fairly standard server computing facilites using a p-value that filters down to 10.000-20.000 variants is reasonable.

```
genotype <-generate.genotype(<input ped>,<input tped>,gwas.id=<selected list of id>,gwas.p=<p-values of input SNPs>)
```

***warning*** If you have more than about 1 million SNPs you must either fast.read = F which will slow down the loading time significantly.  You can also increase your stacklimit using ulimit in the command line,but do this only if you know what you are doing. 


***LD-Filtetring***

There is the option of applying LD filtering using the LD_blocks function after generating the genotype matrix.
Note that this requires that the input data is sorted by chromosome and coordinate to work properly. Read the 
manual discription of the function for more detail.
```
LD_genotype<-LD_blocks(genotype)
genotype <- LD_genotype$genotype
```

***Epistatic Analysis***

After generating the genotypes file it is recomended to run a test run to estimate run time
of the epistatic interaction calculation based on available computing setup:
```
epistatic.correlation(<phenotype dataframe>, genotype,threads = <number of cores available> ,test=T)
```

This will give you an order of magnitude of the expected run time given your input, but not exact time. The next step is to run the analysis:
We recommend using simple=F for better results:
```
correlations<-epistatic.correlation(<phenotype dataframe>, genotype,parallel = <number of cores available> ,test=F,simple=F)
```
Once you have calculated epistatic correlations you can get a coarse grained overview of the results using
the genome.interaction() function:
```
genome.interaction(<input tped file>, correlations)
```

Finally to create modules of interacting markers use the generate.modules function:
```
modules <-generate.modules(correlations)
```

The modules object includes a large range of outputs from the network analysis. 

***warning*** If you have epistatic correlation coefficients values that are strong outliers this can heavily affect
this step or even make it fail. Thus it is recommended to set those coefficients to 0, for example if we only want 
coefficients between 1000 and -1000:
```
correlations$Coefficients[correlations$Coefficients > 1000 |  correlations$Coefficients < -1000] <-0
```
The range of desired values is individual depending on the properties of each dataset.

# Contact
Victor A. O Carmelo, vaocar@bioinformatics.dtu.dk
Quantitative Systems Genetics Group
Department of Bio and Health Informatics
Technical University of Denmark

# References

Lisette J.A. Kogelman and Haja N.Kadarmideen (2014). 
Weighted Interaction SNP Hub (WISH) network method for building genetic 
Networks for complex diseases and traits using whole genome genotype data. 
BMC Systems Biology 8(Suppl 2):S5.  
http://www.biomedcentral.com/1752-0509/8/S2/S5.
