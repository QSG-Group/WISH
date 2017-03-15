partial_correlations_core <- function(genotype_shared,genotype_rev,phenotype,coords,model=1,size){
  n=dim(genotype_shared)[2]
  data_matrix <- matrix(0,nrow = 2*(coords[2]-coords[1]+1),ncol=dim(genotype_shared)[2])
  matrix_row <- 0
  if (model==1){
    for (i in coords[1]:coords[2]){
      matrix_row <- matrix_row+1
      if (i < n){
        for (j in (i+1):n){
          tmp_model = fastLm(phenotype ~ I(genotype_shared[,i])+I(genotype_shared[,j])+I(genotype_shared[,i]*genotype_shared[,j]))
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

partial_correlations_test <- function(phenotype,coords,model=1,size){
  n=size
  data_matrix <- matrix(0,nrow = 2*(coords[2]-coords[1]+1),ncol=size)
  matrix_row <- 0
  if (model==1){
    for (i in coords[1]:coords[2]){
      matrix_row <- matrix_row+1
      if (i < n){
        for (j in (i+1):n){
          tmp_model = fastLm(phenotype ~ I(genotype1[,i])+I(genotype1[,j])+I(genotype1[,i]*genotype1[,j]))
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

epistatic.correlation_test <- function(phenotype,genotype,parallel=1,test=T,simple=T){ 
  phenotype < as.matrix(phenotype)
  n<-ncol(genotype)
  coords<-triangular_split(n,parallel)
  if(is.data.frame(genotype)){
    genotype[] <- lapply(genotype, as.numeric)
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
    x <- as.big.matrix(x = genotype, type = "double", 
                       separated = FALSE, 
                       backingfile = "genotype.file", 
                       descriptorfile = "genotype.desc")
    # get a description of the matrix
    mdesc <- describe(x)
    cl <- makeCluster(parallel)
    registerDoParallel(cl = cl)
    clusterExport(cl,c("partial_correlations_core"))
    message("Running Test")
    message("Estimating run time based on ~100.000 models")
    start.time <- Sys.time()
    test_coords<-triangular_split(1000,parallel)
    snp_matrix <- foreach(j = 1:parallel, .combine='rbind', .verbose=F) %dopar% {
      require(bigmemory)
      require(RcppEigen)
      genotype_shared <- attach.big.matrix("genotype.desc")
      subset <- partial_correlations_core(genotype_shared,1,phenotype,test_coords[j,],model=1)
      return(subset)
    #   n=1000
    #   data_matrix <- matrix(0,nrow = 2*(coords[j,2]-coords[j,1]+1),ncol=1000)
    #   matrix_row <- 0
    #   for (i in coords[j,1]:coords[j,2]){
    #     matrix_row <- matrix_row+1
    #     if (i < n){
    #       for (j in (i+1):n){
    #         tmp_model = fastLm(phenotype ~ I(genotype_shared[,i])+I(genotype_shared[,j])+I(genotype_shared[,i]*genotype_shared[,j]))
    #         data_matrix[(matrix_row*2-1):(matrix_row*2),j]<-c(tmp_model$coefficients[length(tmp_model$coefficients)],summary(tmp_model)$coefficients[dim(summary(tmp_model)$coefficients)[1],4])
    #       }
    #     }
    #   }
    #   return(data_matrix)
    # }
    # message("Running Test")
    # message("Estimating run time based on ~100.000 models")
    # start.time <- Sys.time()
    # genotype_list <- list() 
    # test_coords<-triangular_split(316,parallel)
    #  for (i in 1:dim(test_coords)[1]){
    #   genotype_list[[i]] <- genotype[,test_coords[i,1]:test_coords[i,2]]
    #   print(length(genotype_list))
    #   subset <- partial_correlations_test(genotype_list[[i]],genotype_list[[i]],phenotype,test_coords[i,],316,model=1)
    #   }
    # #snp_matrix <- foreach(j = 1:parallel, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
    #subset <- partial_correlations_test(genotype_list[[j]],genotype_list[[j]],phenotype,dim(genotype_list[[j]])[2],316,model=1)
    #return(subset)
    }
  }
  end.time <- Sys.time()
  time<-as.numeric(end.time-start.time,units="hours")
  model_time<-(((n^2-n)/2)/((316^2-316)/2))*time
  model_time <- round(model_time,digits = 2)
  model_time<-as.character(model_time)
  estimate<-paste(paste("The estimated run time for the simple model is",model_time),"hours",sep=" ")
  message(estimate)
  #snp_matrix <- foreach(j = 1:parallel, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
  # subset <- partial_correlations_test(genotype[,1:316],genotype_rev[,1:316],phenotype,test_coords[j,],model=2)
  #return(subset)
  #}
  #end.time <- Sys.time()
  #time<-as.numeric(end.time-start.time,units="hours")
  #model_time<-(((n^2-n)/2)/((316^2-316)/2))*time
  #model_time <- round(model_time,digits = 2)
  #model_time<-as.character(model_time)
  #estimate<-paste(paste("The estimated run time for the full model is",model_time),"hours",sep=" ")
  #message(estimate)
  epi_cor <- snp_matrix[seq(1,nrow(snp_matrix)-1,2),]
  epi_pvalue <-   snp_matrix[seq(2,nrow(snp_matrix),2),]
  epi_cor_t <- t(epi_cor)
  epi_pvalue_t <- t(epi_pvalue)
  diag(epi_cor_t) <- 1
  diag(epi_pvalue_t) <- 1
  epi_cor_t[upper.tri(epi_cor_t)]<- epi_cor[upper.tri(epi_cor)]
  epi_pvalue_t[upper.tri(epi_pvalue_t)]<- epi_pvalue[upper.tri(epi_pvalue)]
  colnames(epi_pvalue_t) <- colnames(genotype)
  rownames(epi_pvalue_t) <- colnames(genotype)
  colnames(epi_cor_t) <- colnames(genotype)
  rownames(epi_cor_t) <- colnames(genotype)
  output <-list(epi_pvalue_t,epi_cor_t)
  names(output)<-c("Pvalues","Coefficients")
  return(output)
}
result1<-epistatic.correlation_test(genotype = genotype[1:193,1:1000], phenotype = phenotype$W0_W12 ,parallel = 10, test = T)

result1<-epistatic.correlation(genotype = genotype[1:193,1:1000], phenotype = phenotype$W0_W12 ,threads = 10, test = F, simple = T)
hist(result1$Pvalues)
result1$Coefficients


hist(result$Pvalues)


ids_delete <- read.table("/home/victor/Documents/WISH_files/IDs_delete.txt")
phenotype <- read.table("/home/victor/Documents/ADD_files/pheno_for_victor.txt",header = T)
genotype<-generate.genotype(ped = "/home/victor/Documents/ADD_files/INDICES_1_2.ped",tped = "/home/victor/Documents/ADD_files/INDICES1_1_2.tped")
genotype<-genotype[!(rownames(genotype) %in% ids_delete[,1]),]
tped <- fread("/home/victor/Documents/ADD_files/INDICES1_1_2.tped", data.table = F)




#### Developing Blocks ####

correlation_blocks <-function(genotype,threshold=0.9){
  snp_block_matrix <- as.matrix(rep(0,dim(genotype)[2]))
  rownames(snp_block_matrix) <- colnames(genotype)
  n_snp <- 1
  start <- 2
  n_block <- 1
  snps <- dim(genotype)[2]
  snp_block_matrix[n_snp,] <- n_block
  while(start <= snps){
    if (n_snp == snps){
      snp_block_matrix[snps,1] = n_block
    }
    else{
      matches<-sum(genotype[,n_snp]/genotype[,start] == 1,na.rm = T)+sum((genotype[,n_snp]+genotype[,start]) == 3,na.rm = T)+sum(c(c(genotype[,1] == 1.5)+ c(genotype[,3] == 1.5)) == 1 ,na.rm = T)/2
      total <- sum(!is.na(genotype[,n_snp]+genotype[,start]))
      similarity <- matches/total
      if (similarity >= threshold){
        snp_block_matrix[start,1] = n_block
        start <- start+1
      }
      else {
        n_snp <- start
        n_block <- n_block +1
        snp_block_matrix[start,1] = n_block
        start <- start +1
      }
    }
  }
  return(snp_block_matrix)
} 

blocks<-correlation_blocks(genotype,0.9)
tail(blocks)
woot<-sort((table(blocks)),decreasing = T)
woot[1:100]

correlation_blocks_running <-function(genotype,threshold=0.9){
  snp_block_matrix <- as.matrix(rep(0,dim(genotype)[2]))
  rownames(snp_block_matrix) <- colnames(genotype)
  n_snp <- 1
  start <- 2
  n_block <- 1
  snps <- dim(genotype)[2]
  snp_block_matrix[n_snp,] <- n_block
  similarity <- 0
  while(start <= snps){
    if (n_snp == snps){
      snp_block_matrix[snps,1] = n_block
    }
    else{
      matches<-sum(genotype[,n_snp]/genotype[,start] == 1,na.rm = T)+sum((genotype[,n_snp]+genotype[,start]) == 3,na.rm = T)+sum((genotype[,n_snp] == 1.5 + genotype[,start] == 1.5)==1,na.rm = T)
      total <- sum(!is.na(genotype[,n_snp]+genotype[,start]))
      score<-matches+(total-matches)/2
      similarity <- score/total
      if (similarity/(start-n_snp) >= threshold){
        snp_block_matrix[start,1] = n_block
        start <- start+1
      }
      else {
        n_snp <- start
        n_block <- n_block +1
        snp_block_matrix[start,1] = n_block
        start <- start +1
        similarity <- 0
      }
    }
  }
  return(snp_block_matrix)
} 


blocks<-correlation_blocks_running(genotype)
tail(blocks)


correlation_blocks_network <-function(genotype,threshold=0.9){
  snp_block_matrix <- as.matrix(rep(0,dim(genotype)[2]))
  rownames(snp_block_matrix) <- colnames(genotype)
  n_snp <- 1
  start <- 2
  n_block <- 1
  snps <- dim(genotype)[2]
  snp_block_matrix[n_snp,] <- n_block
  while(start <= snps){
    if (n_snp == snps){
      snp_block_matrix[snps,1] = n_block
    }
    else{
      network_coord <- combn(start:n_snp,2)
      similarity <- rep(0,dim(network_coord)[2])
      for (i in 1:dim(network_coord)[2]){
        matches<-sum(genotype[,network_coord[1,i]]/genotype[,network_coord[2,i]] == 1,na.rm = T)+sum((genotype[,network_coord[1,i]]+genotype[,network_coord[2,i]]) == 3,na.rm = T)+sum(c(c(genotype[,network_coord[1,i]] == 1.5)+ c(genotype[,network_coord[2,i]] == 1.5)) == 1 ,na.rm = T)/2
        total <- sum(!is.na(genotype[,network_coord[1,i]]+genotype[,network_coord[2,i]]))
        similarity[i] <- matches/total
      }
      if (sum(similarity)/(dim(network_coord)[2]) >= threshold){
        snp_block_matrix[start,1] = n_block
        start <- start+1
      }
      else {
        n_snp <- start
        n_block <- n_block +1
        print(c(start,n_block,dim(network_coord)[2]))
        snp_block_matrix[start,1] = n_block
        start <- start +1
        similarity <- 0
      }
    }
  }
  return(snp_block_matrix)
} 

blocks<-correlation_blocks_network(genotype,0.95)
tail(blocks)

matches<-sum(tempgen[,1]/tempgen[,2] == 1,na.rm = T)+sum((tempgen[,1]+tempgen[,2]) == 3,na.rm = T)
score<-matches+(sum(!is.na(tempgen[,1]+tempgen[,2]))-matches)/2
sum(c(c(genotype[,1] == 1.5)+ c(genotype[,3] == 1.5)) == 1 ,na.rm = T)

sum((genotype[,1]+genotype[,3]) == 2.5)
sum((genotype[,1]+genotype[,3]) == 3.5)
