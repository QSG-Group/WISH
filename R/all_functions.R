


#' Import genotype data in the correct format for network construction
#' @description For network construction based on both genomic correlations
#' as well as epistatic interactions a genotype matrix has to be
#' created, consisting of one numeric value per SNP, per individual. This function
#' takes Plink output (1,2-coding) to create the genotype matrix which can be used
#' to calculate genomic correlations or epistatic interaction effects 
#' @usage generate.genotype(ped, tped, gwas_id=tped[,2], pvalue=0.05, id.select=ped[,2],gwas_p=NULL,major_freq=0.95)
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
#' @param major_freq Maximum major allele frequency allowed in each variant. 
#' Default value is 0.95. 
#' @return A genotype dataframe and the corresponding vector of passing snps in a vector.
#' The genotype data frame has a row for each individual and a column
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
#' generate.genotype(ped, tped, gwas_id, gwas_p, pvalue, id.select,gwas_p,major_freq)
#' 
#' @export
#' 
#' 
#' 
generate.genotype <- function(ped,tped,gwas_id=tped[,2],pvalue=0.05,id.select=ped[,2],gwas_p=NULL,major_freq=0.95) {
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
    if (length(c(gwas_id))==length(c(tped[,2]))){
      if (length(c(gwas_id))==length(c(tped[,2])))
      {
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
  #Ensuring that we only get variants with enough variation. We remove variants with no minor alleles or/and with a majore allele frequency over 0.95(default)
  passing_snps <- which((colSums((genotype == 2),na.rm = T)*colSums((genotype == 1),na.rm = T)) > 0 & colSums(genotype == 1,na.rm = T) < (dim(genotype)[1]*major_freq))
  genotype <- genotype[,passing_snps]
  snps <- gwas_id[passing_snps]
  return(list(genotype,snps))
}




#' This function calculates the row coordinates for splitting triangular sub
#' matrices of quadratic matrices into approximately equally sized partitions
#' for use in in dividing correlation calculations into equal size for 
#' parallelization    
#' @description Internal function for splitting triangular matrices into
#' approximately equal parts
#' @usage triangular_split(n, split)
#' @param n Row and Column length of the n by n matrix the triangular matrix
#' originates from
#' @param split Number of partitions to split the triangular matrix in
#' @return A matrix of row coordinates used for splitting
#' @examples
#' triangular_split(1000,5)
#' 
#' @export


triangular_split <- function(n,split) {
  if (split == 1){
    boundaries<-matrix(0,nrow=split,ncol=2)
    boundaries[1,1] <- 1
    boundaries[1,2] <- n
  }
  else {
    total_count<-(n*n-n)/2
    splits <- total_count/split
    total <- splits
    row <- c()
    for (i in 1:(split-1)){
      temp_row<-((2*n-1)-sqrt((2*n-1)^2-8*total))/2
      temp_row <- floor(temp_row)
      row <- c(row,temp_row)
      total <- total+splits
    }
    row<-as.vector(c(0,row,n))
    boundaries<-matrix(0,nrow=split,ncol=2)
    for (i in 1:split){
      boundaries[i,1] <- row[i]+1
      boundaries[i,2] <- row[i+1]
    }
  }
  return(boundaries)
}

#' This function calculates the epistatic correlations in a subset of 
#' a matrix space based on coordiantes
#' @description Internal function for calculating epsitatic correlations
#' in sub-matrices
#' @usage partial_correlations(genotype,genotype_rev,phenotype,coords,model)
#' @param genotype Dataframe with the genotype information, resulting from 
#' the function generate.genotype(). Make sure that the dataframe contains the 
#' same individuals as in the phenotype-file, and that those are in the 
#' same order.
#' @param genotype_rev Same as genotpye but with reversed genotype coding
#' @param phenotype Dataframe with the rows correspinding to the individuals
#' in the analysis,and columns for the different measured phenotypes and 
#' fixed/random factors. Phenotypes should be continous variables. 
#' @param coords Matrix of row split coordinates for subseting input space
#' @param model Specification controlling if MM or Mm directed interaction
#' model is used.
#' @return Epsitatic correlations and P-values for the selected set or subset
#' of the data
#' @examples
#' partial_correlations <- function(genotype,genotype_rev,phenotype,coords,model)
#' 
#' @export

partial_correlations <- function(genotype,genotype_rev,phenotype,coords,model=1){
  n=dim(genotype)[2]
  data_matrix <- matrix(0,nrow = 2*(coords[2]-coords[1]+1),ncol=dim(genotype)[2])
  matrix_row <- 0
  if (model==1){
    for (i in coords[1]:coords[2]){
      matrix_row <- matrix_row+1
      if (i < n){
        for (j in (i+1):n){
          tmp_model = fastLm(phenotype ~ I(genotype[,i])+I(genotype[,j])+I(genotype[,i]*genotype[,j]))
          data_matrix[(matrix_row*2-1):(matrix_row*2),j]<-c(tmp_model$coefficients[length(tmp_model$coefficients)],summary(tmp_model)$coefficients[dim(summary(tmp_model)$coefficients)[1],4])
        }
      }
    }
  }
  if (model==2){
    for (i in coords[1]:coords[2]){
      matrix_row <- matrix_row+1
      if (i < n){
        for (j in (i+1):n){
          tmp_model = fastLm(phenotype ~ I(genotype[,i])+I(genotype[,j])+I(genotype[,i]*genotype_rev[,j]))
          data_matrix[(matrix_row*2-1):(matrix_row*2),j]<-c(tmp_model$coefficients[length(tmp_model$coefficients)],summary(tmp_model)$coefficients[dim(summary(tmp_model)$coefficients)[1],4])
        }
      }
    }
  }
  return(data_matrix)
}


#' Calculate the epistatic interaction effect between SNP pairs to construct a 
#' WISH network using a genotype data frame created from genarate.genotype()
#' @import doParallel
#' @import RcppEigen
#' @description A WISH network can be build based on epistatic interaction 
#' effects between SNP pairs. Those interaction effects are calculated using
#' ASReml and can be used directly in the WISH network construction.
#' @usage epistatic.correlation(phenotype, genotype, parallel,test,simple)
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
#' @param test True or False value indicating if a test run is being perform.
#' If True will calculate the expected time it will take for the full analysis
#' based on calculating 100.000 models with the setting chosen
#' @param simple True or false value indicating if only a major/major and
#' minor/minor directed interaction model are tested (simple=T) or if if 
#' interactions on the major/minor minor axis are tested as well, with the 
#' best one of the two being selected (simple=F).
#' @return A list of two matrices. The first matrix gives the epistatic
#' interaction effects between all the SNP-pairs which were in the input 
#' genotype data) and selected with the pvalue from the GWAS results. 
#' The second matrix are the corresponding pvalues of the parameter 
#' estimates of the epistatic interactions. 
#' @references Lisette J.A. Kogelman and Haja N.Kadarmideen (2014). 
#' Weighted Interaction SNP Hub (WISH) network method for building genetic
#' networks for complex diseases and traits using whole genome genotype data.
#' BMC Systems Biology 8(Suppl 2):S5. 
#' http://www.biomedcentral.com/1752-0509/8/S2/S5.
#' @examples
#' epistatic.correlation(phenotype,genotype,parallel,test,simple)
#' 
#' @export

epistatic.correlation <- function(phenotype,genotype,parallel=1,test=T,simple=T){
  registerDoParallel(parallel)
  phenotype < as.matrix(phenotype)
  n<-ncol(genotype)
  coords<-triangular_split(n,parallel)
  if(is.data.frame(genotype)){
    genotype[] <- lapply(genotype, as.numeric)
  }
  else if (is.matrix(genotype)) {
  }
  else {
    stop("genotype not matrix or dataframe")
  }
  if (simple==F || test==T){
    genotype_rev <- genotype
    decide_1<-(genotype_rev==1)
    decide_2<-(genotype_rev==2)
    genotype_rev[decide_1] <- 2
    genotype_rev[decide_2] <- 1
    rm(decide_1)
    rm(decide_2)
    genotype_rev <- as.data.frame(genotype_rev)
  }
  if (test==T && n > 315) {
    message("Running Test")
    message("Estimating run time based on ~100.000 models")
    start.time <- Sys.time()
    test_coords<-triangular_split(316,parallel)
    snp_matrix <- foreach(j = 1:parallel, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations(genotype[,1:316],genotype_rev[,1:316],phenotype,test_coords[j,],model=1)
      return(subset)
    }
    end.time <- Sys.time()
    time<-as.numeric(end.time-start.time,units="hours")
    model_time<-(((n^2-n)/2)/((316^2-316)/2))*time
    model_time <- round(model_time,digits = 2)
    model_time<-as.character(model_time)
    estimate<-paste(paste("The estimated run time for the simple model is",model_time),"hours",sep=" ")
    message(estimate)
    snp_matrix <- foreach(j = 1:parallel, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations(genotype[,1:316],genotype_rev[,1:316],phenotype,test_coords[j,],model=2)
      return(subset)
    }
    end.time <- Sys.time()
    time<-as.numeric(end.time-start.time,units="hours")
    model_time<-(((n^2-n)/2)/((316^2-316)/2))*time
    model_time <- round(model_time,digits = 2)
    model_time<-as.character(model_time)
    estimate<-paste(paste("The estimated run time for the full model is",model_time),"hours",sep=" ")
   message(estimate)
  }
  else if (test==T && n <= 315){
    message("Data size too small for testing, running normal analysis")
    snp_matrix <- foreach(j = 1:parallel, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations(genotype,genotype_rev,phenotype,coords[j,],model=1)
      return(subset)
    }
    # Running opposite minor/major co-linearity model
    snp_matrix_rev <- foreach(j = 1:parallel, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations(genotype,genotype_rev,phenotype,coords[j,],model=2)
      return(subset)
    }
  }
  else if (test==F && simple==F) {
    snp_matrix <- foreach(j = 1:parallel, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations(genotype,genotype_rev,phenotype,coords[j,],model=1)
      return(subset)
    }
    # Running opposite minor/major co-linearity model
    snp_matrix_rev <- foreach(j = 1:parallel, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations(genotype,genotype_rev,phenotype,coords[j,],model=2)
      return(subset)
    }
    # Transposing and filling out the correlation and pvalue matrix
    epi_cor <- snp_matrix[seq(1,nrow(snp_matrix)-1,2),]
    epi_pvalue <-   snp_matrix[seq(2,nrow(snp_matrix),2),]
    rm(snp_matrix)
    epi_cor_t <- t(epi_cor)
    epi_pvalue_t <- t(epi_pvalue)
    diag(epi_cor_t) <- 1
    diag(epi_pvalue_t) <- 1
    epi_cor_t[upper.tri(epi_cor_t)]<- epi_cor[upper.tri(epi_cor)]
    epi_pvalue_t[upper.tri(epi_pvalue_t)]<- epi_pvalue[upper.tri(epi_pvalue)]
    epi_cor_t_1 <- epi_cor_t
    epi_pvalue_t_1 <- epi_pvalue_t
    epi_cor_t_1[is.na(epi_cor_t_1)] <- 0
    epi_pvalue_t_1[is.na(epi_pvalue_t_1)] <- 1
    # Opposite assumption correlation and pvalue matrix
    epi_cor <- snp_matrix_rev[seq(1,nrow(snp_matrix_rev)-1,2),]
    epi_pvalue <-   snp_matrix_rev[seq(2,nrow(snp_matrix_rev),2),]
    rm(snp_matrix_rev)
    epi_cor_t <- t(epi_cor)
    epi_pvalue_t <- t(epi_pvalue)
    diag(epi_cor_t) <- 1
    diag(epi_pvalue_t) <- 1
    epi_cor_t[upper.tri(epi_cor_t)]<- epi_cor[upper.tri(epi_cor)]
    epi_pvalue_t[upper.tri(epi_pvalue_t)]<- epi_pvalue[upper.tri(epi_pvalue)]
    epi_cor_t_2 <- epi_cor_t
    epi_pvalue_t_2 <- epi_pvalue_t
    epi_cor_t_2[is.na(epi_cor_t_2)] <- 0
    epi_pvalue_t_2[is.na(epi_pvalue_t_2)] <- 1
    #Picking the correct model assumption based on lowest pvalue
    decider_matrix <- epi_pvalue_t_1-epi_pvalue_t_2
    epi_pvalue_t <- epi_pvalue_t_1
    epi_pvalue_t[0 < decider_matrix] <- epi_pvalue_t_2[0 < decider_matrix]
    epi_cor_t <- epi_cor_t_1
    epi_cor_t[0 < decider_matrix] <- epi_cor_t_2[0 < decider_matrix]
    return(list(epi_pvalue_t,epi_cor_t))
  }
  else if (test==F && simple==T) {
    snp_matrix <- foreach(j = 1:parallel, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations(genotype,genotype_rev,phenotype,coords[j,],model=1)
      return(subset)
    }
    epi_cor <- snp_matrix[seq(1,nrow(snp_matrix)-1,2),]
    epi_pvalue <-   snp_matrix[seq(2,nrow(snp_matrix),2),]
    epi_cor_t <- t(epi_cor)
    epi_pvalue_t <- t(epi_pvalue)
    diag(epi_cor_t) <- 1
    diag(epi_pvalue_t) <- 1
    epi_cor_t[upper.tri(epi_cor_t)]<- epi_cor[upper.tri(epi_cor)]
    epi_pvalue_t[upper.tri(epi_pvalue_t)]<- epi_pvalue[upper.tri(epi_pvalue)]
    return(list(epi_pvalue_t,epi_cor_t))
  }
}


