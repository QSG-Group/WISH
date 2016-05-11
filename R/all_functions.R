#library("GenABEL")
#library("asreml")
#library("doParallel")
#library("sna")

#' Import genotype data in the correct format for network construction
#' @export
#' @description For network construction based on both
#' @usage data.import(ped, tped, gwas_id, gwas_p, pvalue, id.select)
#' genomic correlations aswell epistatic interactions a genotype matrix has to be
#' created, consisting of one numeric value per SNP, per individual. This function
#' takes Plink output (1,2-coding) to create the genotype matrix which can be used
#' to calculate genomic correlations or epistatic interaction effects 
#' @param ped The ped file (.ped) is an input file from Plink: The PED file is a
#' white-space (space or tab) delimited file: the first six columns are mandatory:
#' Family ID, Idividual ID, Paternal ID, Maternal ID, 
#' Sex (1=male; 2=female;other=unknown) and Phenotype. The IDs are alphanumeric: 
#' the combination of family and individual ID should uniquely identify a person.
#' A PED file must have1 and only 1 phenotype in the sixth column.
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
#' @param gwas_p A vector of the p-values corresponding to the gwas_id vector
#' @param pvalue A value for the cutoff of the SNPs which should be remained 
#' in the matrix, based on the pvalue resulting from the GWAS
#' @param id.select If requested, asubset of individuals can be 
#' selected (e.g. extremes). If nothing inserted, all individuals are in the
#' output
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
#' data.import(ped, tped, gwas_id, gwas_p, pvalue, id.select)
#' 
#'
data.import <- function(ped, tped,gwas_id,gwas_p,pvalue=0.05,id.select=ped[,2]){
  a.s <- ped[,c(2,7:ncol(ped))]
  snp <- rep((as.character(tped[,2])),each=2)
  snp <- append(snp, "ID", after=0)
  colnames(a.s) <- c(snp)
  gwas.select <- cbind(gwas_id,gwas_p)
  gwas.select <- as.data.frame(gwas.select)
  gwas.select$gwas_p <- as.numeric(as.character(gwas.select$gwas_p))
  p.select <- gwas.select[gwas.select$gwas_p < pvalue,]
  snp.select <- p.select$gwas_id
  snps<- rep(snp.select, each=2)
  snps <- append(as.character(snps), "ID", after=0)
  a.snps <- a.s[,snps]
  rownames(a.snps) <- a.snps[,1]
  input <- a.snps[,2:ncol(a.snps)]
  ind <- data.frame(matrix(c((1:ncol(input)),rep(NA, 2-ncol(input)%%2)),byrow=F,nrow=2))
  nonna <- ind[,sapply(ind, function(x) all(!is.na(x)))]
  genotype <- do.call(cbind,lapply(nonna,function(i)rowMeans(input[,i])))
  x <- colnames(input[1:(ncol(input)-1)])
  colnames(genotype) <- x[c(T,F)]
  genotype <- replace(genotype, genotype==0,"NA")
  return(as.data.frame(genotype))
}

#' Calculate the epistatic interaction effect between SNP pairs to construct a 
#' WISH network
#' @export
#' @description A WISH network can be build based on epistatic interaction 
#' effects between SNP pairs. Those interaction effects are calculated using
#' ASReml and can be used directly in the WISH network construction.
#' @usage epistatic.correlation(gwas, pvalue, phenotype, genotype, parallel)
#' ASReml-R: library("asreml")
#' package doParralel: library("doParallel")
#' package sna: library("sna")
#' @import asreml
#' @import doParallel
#' @import sna
#' @param gwas_id A vector of all SNPs in the GWAS
#' @param gwas_p A vector of the p-values corresponding to the gwas_id vector
#' @param pvalue A value for the cutoff of the SNPs which should be remained 
#' in the matrix, based on the pvalue resulting from the GWAS
#' @param phenotype Dataframe with on the rows the individuals in the analysis,
#' and columns for the different measured phenotypes and fixed/random factors
#' (e.g. sex)
#' @param genotype Dataframe with the genotype information, resulting from 
#' the function data.import. Make sure that the dataframe contains the same
#' individuals as in the phenotype-file, and that those are in the same order.
#' @param parallel Number of cores to use for parallel execution in the function 
#' registerDoParallel()
#' @return The resulting matrix gives the epistatic interaction effects between
#' all the SNP-pairs which were in the input (genotype data) and selected with
#' the pvalue from the GWAS results. 
#' @references Lisette J.A. Kogelman and Haja N.Kadarmideen (2014). 
#' Weighted Interaction SNP Hub (WISH) network method for building genetic
#' networks for complex diseases and traits using whole genome genotype data.
#' BMC Systems Biology 8(Suppl 2):S5. 
#' http://www.biomedcentral.com/1752-0509/8/S2/S5.
#' @examples
#' epistatic.correlation(gwas, pvalue, phenotype, genotype, parallel)
#' 
#' 
epistatic.correlation <- function(gwas_id,gwas_p,pvalue=10E-5,phenotype,genotype,parallel=20 ){
  registerDoParallel(parallel)
  gwas <- cbind(gwas_id,gwas_p)
  gwas <- as.data.frame(gwas)
  gwas$gwas_p <- as.numeric(as.character(gwas$gwas_p))
  gwas <- gwas[gwas$gwas_id%in%colnames(genotype),]
  cand <- gwas[gwas$gwas_p < pvalue,]
  snp_matrix <- matrix(NA, nrow=nrow(cand),ncol=nrow(cand))
  rownames(snp_matrix)<- cand$gwas_id
  colnames(snp_matrix)<- cand$gwas_id
  cand <- which(gwas$gwas_p<= pvalue)  
  for (i in 1:(nrow(snp_matrix)-1)) {
    snp_matrix[i, (i+1):nrow(snp_matrix)] <- foreach(j = (i+1):nrow(snp_matrix), .combine='rbind', .inorder=T, .verbose=F) %dopar% { 
      phenotype$SNP1=as.numeric((genotype[,cand[i]]));
      phenotype$SNP2=as.numeric((genotype[,cand[j]]));
      tmp_correlation = asreml(OI~1+SNP1+SNP2+SNP1*SNP2,na.method.X="omit", na.method.Y="omit",data=phenotype) ## <- this function needs to be adjusted according to your trait: trait=trait of interest, fixed and random effects can be added (e.g. sex)
      return(tmp_correlation$coefficients$fixed[1])
    } 
  }
  new <- t(snp_matrix)
  diag(new) <- 1 
  new[upper.tri(new)]<- snp_matrix[upper.tri(snp_matrix)]
  return(new)
}  
