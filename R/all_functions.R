

#Loading Packages

#' Import genotype data in the correct format for network construction
#' @description For network construction based on both #' genomic correlations
#' as well as epistatic interactions a genotype matrix has to be
#' created, consisting of one numeric value per SNP, per individual. This function
#' takes Plink output (1,2-coding) to create the genotype matrix which can be used
#' to calculate genomic correlations or epistatic interaction effects 
#' @usage generate.genotype(ped, tped, gwas_id, pvalue=0.05, id.select=ped[,2],gwas_p=NULL)
#' @param ped The ped file (.ped) is an input file from Plink: The PED file is a
#' white-space (space or tab) delimited file: the first six columns are mandatory:
#' Family ID, Idividual ID, Paternal ID, Maternal ID, 
#' Sex (1=male; 2=female;other=unknown) and Phenotype. The IDs are alphanumeric: 
#' the combination of family and individual ID should uniquely identify a person.
#' A PED file must have 1 and only 1 phenotype in the sixth column.
#' The phenotype can be either a quantitative trait or an affection status 
#' column: PLINK will automatically detect which type
#' (i.e. based on whether a value other than 0, 1, 2 or the missing genotype 
#' code is observed). SNPs are 1,2-coded (1 for major allele,2 for minor allele) 
#' For more information: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
#' @param tped The tped file (.tped) is a transposed ped file, from Plink. 
#' This file contains the SNP and genotype information where one row is a SNP.
#' The first 4 columns of a TPED file are the same as a 4-column MAP file.
#' Then all genotypes are listed for all individuals for each particular SNP on 
#' each line. Again, SNPs are 1,2-coded.
#' @param gwas_id  A vector of all SNPs in the GWAS
#' @param pvalue A value for the cutoff of the SNPs which should be remained 
#' in the matrix, based on the pvalue resulting from the GWAS. Default value
#' is 0.05
#' @param id.select If requested, a subset of individuals can be 
#' selected (e.g. extremes). If nothing inserted, all individuals are in the
#' output
#' @param gwas_p **optional** A vector of the p-values corresponding to 
#' the gwas_id vector. If assigned, will select snps based on the pvalue
#' parameter with a default value of 0.05.
#' @return A genotype dataframe with a row for each individual and a column
#'  for each SNP. SNPs are 1,1.5,2 coded: 1 for homozygous for the major 
#'  allele, 1.5 for heterozygous, and 2 for homozygous for the minor allele. 
#'  Missing values are NA coded. 
#' @details There is so much to be said
#' @references Lisette J.A. Kogelman and Haja N.Kadarmideen (2014). 
#' Weighted Interaction SNP Hub (WISH) network method for building genetic
#' networks for complex diseases and traits using whole genome genotype data.
#' BMC Systems Biology 8(Suppl 2):S5. 
#' http://www.biomedcentral.com/1752-0509/8/S2/S5.
#' @examples
#' generate.genotype(ped, tped, gwas_id, gwas_p, pvalue, id.select)
#' 
#' @export
generate.genotype <- function(ped,tped,gwas_id=tped[,2],pvalue=0.05,id.select=ped[,2],gwas_p=NULL) {
  if(is.null(gwas_p)){
    genotype <- matrix(nrow=length(c(id.select)),ncol=length(c(gwas_id)))
    rownames(genotype) <- id.select
    colnames(genotype) <- gwas_id
    if (length(c(gwas_id))==length(c(tped[,2]))){
      snps <- c(1:dim(tped)[1])
    }
    else {
      snps<-which(tped[,2]%in%gwas_id)  
    }
    if (length(c(id.select))==length(c(ped[,2]))){
      ids <- c(1:dim(ped)[1])
      ped_trim <- as.matrix(ped[ids,c(rep(2*snps,each=2)-(1:(2*length(snps)))%%2+6)])
    }
    else {
      ids<-which(ped[,2]%in%id.select)
      ped_trim <- as.matrix(ped[ids,c(sort(rep(2*snps,each=2)-(1:(2*length(snps)))%%2)+6)])
    }
    ped_trim[ped_trim==0] <- NA
    for (i in 1:(dim(genotype)[2])){
      genotype[,i] <- rowMeans((ped_trim[,c(2*i-1,2*i)]))
    }
  }
  if(!(is.null(gwas_p))){
    gwas_id <- as.vector(gwas_id)
    gwas_id <- gwas_id[as.vector(gwas_p) <= pvalue]
    genotype <- matrix(nrow=length(c(id.select)),ncol=length(c(gwas_id)))
    rownames(genotype) <- id.select
    colnames(genotype) <- gwas_id
    if (length(c(gwas_id))==length(c(tped[,2]))){if (length(c(gwas_id))==length(c(tped[,2]))){
      snps <- c(1:dim(tped)[1])
    }
      else {
        snps<-which(tped[,2]%in%gwas_id)  
      }
      if (length(c(id.select))==length(c(ped[,2]))){
        ids <- c(1:dim(ped)[1])
        ped_trim <- ped[ids,c(rep(2*snps,each=2)-(1:(2*length(snps)))%%2+6)]
      }
      else {
        ids<-which(ped[,2]%in%id.select) 
        ped_trim <- ped[ids,c(sort(rep(2*snps,each=2)-(1:(2*length(snps)))%%2)+6)]
      }
      snps <- c(1:dim(tped)[1])
    }
    else {
      snps<-which(tped[,2]%in%gwas_id)  
    }
    if (length(c(id.select))==length(c(ped[,2]))){
      ids <- c(1:dim(ped)[1])
      ped_trim <- as.matrix(ped[ids,c(rep(2*snps,each=2)-(1:(2*length(snps)))%%2+6)])
    }
    else {
      ids<-which(ped[,2]%in%id.select) 
      ped_trim <- as.matrix(ped[ids,c(sort(rep(2*snps,each=2)-(1:(2*length(snps)))%%2)+6)])
    }
    ped_trim[ped_trim==0] <- NA
    for (i in 1:(dim(genotype)[2])){
      genotype[,i] <- rowMeans((ped_trim[,c(2*i-1,2*i)]))
    }
  }
  return(genotype)
}

#' Calculate the epistatic interaction effect between SNP pairs to construct a 
#' WISH network using a genotype data frame created from genarate.genotype()
#' @import doParallel
#' @import RcppEigen
#' @description A WISH network can be build based on epistatic interaction 
#' effects between SNP pairs. Those interaction effects are calculated using
#' ASReml and can be used directly in the WISH network construction.
#' @usage epistatic.correlation(phenotype, genotype, parallel)
#' @param phenotype Dataframe with the rows correspinding to the individuals
#' in the analysis,and columns for the different measured phenotypes and 
#' fixed/random factors. Only give one phenotype column at a time. Phenotypes
#' should be continous variables. Make sure that the dataframe contains the same
#' individuals as in the genotype-file, and that those are in the same order.
#' @param genotype Dataframe with the genotype information, resulting from 
#' the function generate.genotype(). Make sure that the dataframe contains the 
#' same individuals as in the phenotype-file, and that those are in the 
#' same order.
#' @param parallel Number of cores to use for parallel execution in the function 
#' registerDoParallel()
#' @return A list of two matrices. The first matrix gives the epistatic
#'  interaction effects between
#' all the SNP-pairs which were in the input (genotype data) and selected with
#' the pvalue from the GWAS results. The second matrix are the corresponding
#' pvalues of the parameter estimates of the epistatic interactions. 
#' @references Lisette J.A. Kogelman and Haja N.Kadarmideen (2014). 
#' Weighted Interaction SNP Hub (WISH) network method for building genetic
#' networks for complex diseases and traits using whole genome genotype data.
#' BMC Systems Biology 8(Suppl 2):S5. 
#' http://www.biomedcentral.com/1752-0509/8/S2/S5.
#' @examples
#' epistatic.correlation(phenotype,genotype,parallel)
#' 
#' @export
epistatic.correlation <- function(phenotype,genotype,parallel=1 ){
  registerDoParallel(parallel)
  snp_matrix <- matrix(NA, nrow=2*(ncol(genotype)),ncol=ncol(genotype))
  rownames(snp_matrix)<- rep(colnames(genotype),each=2)
  colnames(snp_matrix)<- colnames(genotype)
  if(is.data.frame(genotype)){
    genotypen <- genotype
    genotypen[] <- lapply(genotypen, as.numeric)
    genotypen<-as.matrix(genotypen)
  }
  else {
    genotypen <- genotype
  }
  genotype <- as.data.frame(genotype)
  genotype[] <- lapply(genotype, as.factor)
  #The structure here ensures we only calculate one half of the matrix, and return in such a way to get both epistatic interactions and pvalues i
  # one go.
  for (i in 1:(ncol(snp_matrix)-1)) {
    snp_matrix[c(i+(i-1),2*i),(i+1):ncol(snp_matrix)] <- foreach(j = (i+1):ncol(snp_matrix), .combine='cbind', .inorder=T, .verbose=F) %dopar% {
      tmp_model = fastLm(phenotype[,1] ~ genotype[,i]+genotype[,j]+I(genotypen[,i]*genotypen[,j]))
      return(c(tmp_model$coefficients[length(tmp_model$coefficients)],summary(tmp_model)$coefficients[dim(summary(tmp_model)$coefficients)[1],4]))
    }
  }
  # Transposing and filling out the correlation and pvalue matrix
  epi_cor <- snp_matrix[seq(1,nrow(snp_matrix)-1,2),]
  epi_pvalue <-   snp_matrix[seq(2,nrow(snp_matrix),2),]
  epi_cor_t <- t(epi_cor)
  epi_pvalue_t <- t(epi_pvalue)
  diag(epi_cor_t) <- 1
  diag(epi_pvalue_t) <- 1
  epi_cor_t[upper.tri(epi_cor_t)]<- epi_cor[upper.tri(epi_cor)]
  epi_pvalue_t[upper.tri(epi_pvalue_t)]<- epi_pvalue[upper.tri(epi_pvalue)]
  return(list(epi_cor_t,epi_pvalue_t))
}  
