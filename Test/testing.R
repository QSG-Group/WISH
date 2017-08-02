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
genotype<-generate.genotype(ped = "/home/victor/Documents/ADD_files/INDICES_1_2.ped",tped = "/home/victor/Documents/ADD_files/INDICES1_1_2.tped",major.freq = 0.95)
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
  block_coords <- list()
  while(start <= snps){
    # if (n_snp == snps){
    #   snp_block_matrix[snps,1] = n_block
    # }
    #else{
      #matches<-sum(genotype[,n_snp]/genotype[,start] == 1,na.rm = T)+(sum(c(c(genotype[,n_snp] == 1.5)+ c(genotype[,start] == 1.5)) == 1 ,na.rm = T)/2)
      total <- sum(!is.na(genotype[,n_snp]+genotype[,start]))
      #similarity <- matches/total
      p_AB <- (sum(genotype[,n_snp]*genotype[,start] == 4,na.rm = T)+(sum(genotype[,n_snp]+ genotype[,start] == 3.5 ,na.rm = T)/2)+(sum(genotype[,n_snp]+ genotype[,start] == 3 ,na.rm = T)/2))/total
      p_A <- (sum(genotype[,n_snp] == 2,na.rm = T)+ sum(genotype[,n_snp] == 1.5,na.rm = T)/2)/total
      p_B <- (sum(genotype[,start] == 2,na.rm = T)+ sum(genotype[,start] == 1.5,na.rm = T)/2)/total
      #print(c(p_A,p_B,p_AB,start))
      D <- p_AB-p_A*p_B
      #print(D)
      #corr <- D^2/(p_A*p_B*(1-p_A)*(1-p_B))
      if (D < 0){
        Dd <- max(-p_AB,-(1-p_A)*(1-p_B))
      }
      else {
        Dd <-min(p_A*(1-p_B),p_B*(1-p_A))
      }
      similarity<-D/Dd
      #print(similarity)
      if (similarity >= threshold){ 
        snp_block_matrix[start,1] = n_block
        start <- start+1
        if (start > snps){
          block_coords[[n_block]] <- c(n_snp,(start-1))
        }
      } 
      else {
        if (start == snps){
          block_coords[[n_block]] <- c(n_snp,(start-1))
          n_block <- n_block + 1
          snp_block_matrix[snps,1] <- n_block
          block_coords[[n_block]] <- c(start,start)
          start <- start + 1
        }
        else{
        block_coords[[n_block]] <- c(n_snp,(start-1))
        n_snp <- start
        n_block <- n_block+1
        snp_block_matrix[start,1] = n_block
        start <- start +1
        if (start%%1000 == 0){
          print(start)
        }
        }
      # if (start == snps){
      # 
      # }
      }
    #}
  }
  genotype <- as.matrix(genotype)
  new_genotype <- matrix(0L, nrow = dim(genotype)[1], ncol = n_block)
  counter <- 0
  message("Creating Blocks")
  for (coords in block_coords){
    counter <- counter + 1
    # if (counter%%1000 == 0 ){
    #   print(counter)
    # }
    if (coords[1] == coords[2]){
      new_genotype[,counter] <- genotype[,coords[1]]
    }
    else {
    new_genotype[,counter] <- rowMeans(genotype[,coords[1]:coords[2]],na.rm = T) 
    }
  }
  output<-list(new_genotype,block_coords,snp_block_matrix)
  names(output)<-c("genotype","block_coords","snp_id_blocks")
  return(output)
}

genotype_cut<-correlation_blocks(genotype)
genotype_cut$genotype
blocks1<-genotype_cut[[1]]
blocks2<-genotype_cut[[2]]
blocks3<-genotype_cut[[3]]

head(blocks2)

head(blocks3)


<dim(genotype_cut[[1]])



correlation_blocks2 <-function(genotype,threshold=0.9){ 
  snp_block_matrix <- as.matrix(rep(0,dim(genotype)[2]))
  rownames(snp_block_matrix) <- colnames(genotype)
  n_snp <- 1
  start <- 2
  n_block <- 1
  snps <- dim(genotype)[2]
  snp_block_matrix[n_snp,] <- n_block
  block_coords <- list()
  while(start <= snps){
    if (n_snp == snps){
      snp_block_matrix[snps,1] = n_block
    }
    else{
      #matches<-sum(genotype[,n_snp]/genotype[,start] == 1,na.rm = T)+(sum(c(c(genotype[,n_snp] == 1.5)+ c(genotype[,start] == 1.5)) == 1 ,na.rm = T)/2)
      total <- sum(!is.na(genotype[,n_snp]+genotype[,start]))
      #similarity <- matches/total
      p_MM <- 2*(sum(genotype[,n_snp]*genotype[,start] == 4,na.rm = T)+(sum(genotype[,n_snp]+ genotype[,start] == 3 ,na.rm = T)))
      p_M <- sum(genotype[,n_snp]+genotype[,start] == 3.5 ,na.rm = T)
      similarity <- (p_MM+p_M)/(total*2)
      print(similarity)
      if (similarity >= threshold){
        snp_block_matrix[start,1] = n_block
        start <- start+1
      }
      else {
        block_coords[[n_block]] <- c(n_snp,start)
        n_snp <- start
        n_block <- n_block +1
        snp_block_matrix[start,1] = n_block
        start <- start +1
        if (start%%1000 == 0){
          print(start)
        }
      }
    } 
  }
  genotype <- as.matrix(genotype)
  new_genotype <- matrix(0L, nrow = dim(genotype)[1], ncol = n_block)
  counter <- 0
  message("Creating Blocks")
  for (coords in block_coords){
    counter <- counter + 1
    if (counter%%1000 == 0 ){
      print(counter)
    }
    if (coords[1] == coords[2]){
      new_genotype[,counter] <- genotype[,coords[1]]
    }
    else {
      new_genotype[,counter] <- rowMeans(genotype[,coords[1]:coords[2]],na.rm = T) 
    }
  }
  output<-list(new_genotype,block_coords)
  names(output)<-c("genotype","block_coords","snp_id_blocks")
  return(list(output)
}

genotype_cut<-correlation_blocks(genotype)
genotype_cut[[2]]

genotype_cut[,1392]

dim(genotype_cut)


system.time(blocks<-correlation_blocks(genotype,0.9))
dim(blocks)
tail(blocks)
genotype[,blocks[,1] == 1]
sum(blocks[,1] == 1)
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
      matches<-sum(genotype[,n_snp]/genotype[,start] == 1,na.rm = T) + sum((genotype[,n_snp]+genotype[,start]) == 3,na.rm = T)+sum(((genotype[,n_snp] == 1.5 + genotype[,start] == 1.5)==1,na.rm = T)
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


correlation_blocks_network <-function(genotype,threshold=0.9,max_block_size=1000){
  snp_block_matrix <- as.matrix(rep(0,dim(genotype)[2]))
  rownames(snp_block_matrix) <- colnames(genotype)
  n_snp <- 1
  start <- 2
  n_block <- 1
  similarity <- 0
  network_size <- 0
  snps <- dim(genotype)[2]
  snp_block_matrix[n_snp,] <- n_block
  while(start <= snps){
    if (n_snp == snps){
      snp_block_matrix[snps,1] = n_block
    }
    else{
      temp_sim <- 0
      for (i in n_snp:(start-1)){
        matches1<-sum(genotype[,start]/genotype[,i] == 1,na.rm = T)+sum(c(c(genotype[,start] == 1.5)+ c(genotype[,i] == 1.5)) == 1 ,na.rm = T)/2
        matches2<-sum(genotype[,start]*genotype[,i] == 2,na.rm = T)+sum(c(c(genotype[,start] == 1.5)+ c(genotype[,i] == 1.5)) == 1 ,na.rm = T)/2
        matches <- max(matches1,matches2)
        total <- sum(!is.na(genotype[,start]+genotype[,i]))
        p_AB <- (sum(genotype[,n_snp]*genotype[,start] == 4,na.rm = T)+(sum(genotype[,n_snp]+ genotype[,start] == 3.5 ,na.rm = T)/2)+(sum(genotype[,n_snp]+ genotype[,start] == 3 ,na.rm = T)/2))/total
        p_A <- (sum(genotype[,n_snp] == 2,na.rm = T)+ sum(genotype[,n_snp] == 1.5,na.rm = T)/2)/total
        p_B <- (sum(genotype[,start] == 2,na.rm = T)+ sum(genotype[,start] == 1.5,na.rm = T)/2)/total
        #print(c(p_A,p_B,p_AB,start))
        D <- p_AB-p_A*p_B
        #print(D)
        #corr <- D^2/(p_A*p_B*(1-p_A)*(1-p_B))
        if (D < 0){
          Dd <- max(-p_AB,-(1-p_A)*(1-p_B))
        }
        else {
          Dd <-min(p_A*(1-p_B),p_B*(1-p_A))
        }
        temp_sim<-D/Dd
        network_size <- network_size + 1
        temp_sim <- temp_sim + matches/total
      }
      similarity <- similarity + temp_sim
      if (similarity/network_size >= threshold && network_size < 1001){
        snp_block_matrix[start,1] = n_block
        start <- start+1
      }
      else {
        n_snp <- start
        n_block <- n_block +1
        snp_block_matrix[start,1] = n_block
        start <- start +1
        if (start %% 1000 == 0){
        print(c(start,network_size,n_block))
        }
        similarity <- 0
        network_size <- 0
      }
    }
  }
  return(snp_block_matrix)
} 

blocks<-correlation_blocks_network(genotype,0.9)
tail(blocks)


matches<-sum(genotype[,1001]/genotype[,1000] == 1,na.rm = T)+sum(c(c(genotype[,1001] == 1.5)+ c(genotype[,1000] == 1.5)) == 1 ,na.rm = T)/2

matches<-sum(tempgen[,1]/tempgen[,2] == 1,na.rm = T)+sum((tempgen[,1]+tempgen[,2]) == 3,na.rm = T)
score<-matches+(sum(!is.na(tempgen[,1]+tempgen[,2]))-matches)/2
sum(c(c(genotype[,1] == 1.5)+ c(genotype[,3] == 1.5)) == 1 ,na.rm = T)

sum((genotype[,1]+genotype[,3]) == 2.5)
sum((genotype[,1]+genotype[,3]) == 3.5)


### New matrix Split ###

#' Calculate the epistatic interaction effect between SNP pairs to construct a 
#' WISH network using a genotype data frame created from genarate.genotype()
#' @import doParallel
#' @import foreach
#' @import RcppEigen
#' @description A WISH network can be built based on epistatic interaction 
#' effects between SNP pairs. Those interaction effects are calculated using
#' linear models. 
#' @usage epistatic.correlation(phenotype, genotype, parallel,test,simple)
#' @param phenotype Dataframe with the rows correspinding to the individuals
#' in the analysis,and columns for the different measured phenotypes and 
#' fixed/random factors. Only give one phenotype column at a time. Phenotypes
#' should be non-categorical continous or discrete/semi-discrete variables. 
#' Make sure that the dataframe contains the same
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

epistatic.correlation_test <- function(phenotype,genotype,threads=1,test=T,simple=T){ 
  registerDoParallel(threads)
  if ( simple == T){ 
    model <- 1
  }
  else {
    model <- 2
  }
    phenotype < as.matrix(phenotype)
  n<-ncol(genotype)
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
    message("Running Test")
    message("Estimating run time based on ~100.000 models")
    start.time <- Sys.time()
    coord_splits <-triangular_split(316,threads)
    snp_matrix <- foreach(j = 1:threads, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations_triangular(genotype[,1:316],genotype_rev[,1:316],phenotype=phenotype,coord_splits[j,],model=model) 
      return(subset)
    }
    end.time <- Sys.time()
    time<-as.numeric(end.time-start.time,units="hours")
    model_time<-(((n^2-n)/2)/((316^2-316)/2))*time
    model_time <- round(model_time,digits = 2)
    model_time<-as.character(model_time)
    estimate<-paste(paste("The estimated run time for the simple model is",model_time),"hours",sep=" ")
    message(estimate)
    coord_splits <-triangular_split(316,threads)
    snp_matrix <- foreach(j = 1:threads, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations_triangular(genotype[,1:316],genotype_rev[,1:316],phenotype=phenotype,coord_splits[j,],model=2) 
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
    coord_splits <-triangular_split(dim(genotype)[2],threads)
    snp_matrix <- foreach(j = 1:threads, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations_triangular(genotype,genotype_rev,phenotype=phenotype,coord_splits[j,],model=model) 
      return(subset)
    }
    # Running opposite minor/major co-linearity model
    coord_splits <-triangular_split(dim(genotype)[2],threads)
    snp_matrix_rev <- foreach(j = 1:threads, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations_triangular(genotype,genotype_rev,phenotype=phenotype,coord_splits[j,],model=model) 
      return(subset)
    }
  }
  else if (test==F && simple==F) {
    coord_splits <-triangular_split(dim(genotype)[2],threads)
    snp_matrix <- foreach(j = 1:threads, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations_triangular(genotype,genotype_rev,phenotype=phenotype,coord_splits[j,],model=model) 
      return(subset)
    }
    # Running opposite minor/major co-linearity model
    coord_splits <-triangular_split(dim(genotype)[2],threads)
    snp_matrix_rev <- foreach(j = 1:threads, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations_triangular(genotype,genotype_rev,phenotype=phenotype,coord_splits[j,],model=model) 
      return(subset)
    }
  }
  if (test == F && simple==F || (test==T && n <= 315)){
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
    colnames(epi_pvalue_t) <- colnames(genotype)
    rownames(epi_pvalue_t) <- colnames(genotype)
    colnames(epi_cor_t) <- colnames(genotype)
    rownames(epi_cor_t) <- colnames(genotype)
    output <-list(epi_pvalue_t,epi_cor_t)
    names(output)<-c("Pvalues","Coefficients")
    return(output)
  }
  else if (test==F && simple==T){ 
    coord_splits <-triangular_split(dim(genotype)[2],threads)
    snp_matrix <- foreach(j = 1:threads, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations_triangular(genotype,genotype_rev,phenotype=phenotype,coord_splits[j,],model=model) 
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
    colnames(epi_pvalue_t) <- colnames(genotype)
    rownames(epi_pvalue_t) <- colnames(genotype)
    colnames(epi_cor_t) <- colnames(genotype)
    rownames(epi_cor_t) <- colnames(genotype)
    output <-list(epi_pvalue_t,epi_cor_t)
    names(output)<-c("Pvalues","Coefficients")
    return(output)
  }
} 

system2("ulimit","-s 16384")
system("ulimit -s")
ids_delete <- read.table("/data/nbg153/ADHD/WISH_files_newest/IDs_delete.txt")
phenotype <- read.table("//data/nbg153/ADHD/WISH_files_newest/pheno_for_victor.txt",header = T)
genotype<-generate.genotype(ped = "/data/nbg153/ADHD/WISH_files_newest/trimmed_ped",tped = "/data/nbg153/ADHD/WISH_files_newest/trimmed_tped",major.freq = 0.95)
genotype<-genotype[!(rownames(genotype) %in% ids_delete[,1]),]
tped <- fread("/data/nbg153/ADHD/WISH_files_newest/INDICES1_1_2.tped", data.table = F)

new_geno <- cbind(genotype,genotype,genotype)

dim(genotype)

result<-epistatic.correlation_test(phenotype[,3],genotype,threads= 40,test=F,simple=F)



phenotype
triangular_split_test <- function(n,split) {
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

# Creates coordinates splits for rectangular matrices
square_split <- function(threads,rows) { 
  split <- floor(rows/threads)
  coords <- matrix(0,nrow=2,ncol=threads)
  coords <- c()
  for (i in 1:threads){ 
    if (i != threads){
      coords <-c(coords,(i)*split)
    } 
  }
  coords <- c(0,coords,rows) 
  boundaries<-matrix(0,nrow=threads,ncol=2)
  for (i in 1:threads){
    boundaries[i,1] <- coords[i]+1
    boundaries[i,2] <- coords[i+1]
  }
  return(boundaries)
}

square_split(10,1000)

memory_subspace <- function(n,threads,genotype,genotype_rev=0,phenotype,model) { 
  # Split correlation into n^2/2 section for memory usage saving
  registerDoParallel(threads)
  if (n == 1){ 
    coord_splits <-triangular_split_test(dim(genotype)[2],threads)
    snp_matrix <- foreach(j = 1:threads, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations_triangular(genotype,genotype_rev,phenotype=phenotype,coord_splits[j,],model=model) 
      return(subset)
    }
  }
  else {
  square_size <- dim(genotype)[2]
  snp_matrix <- matrix(0,nrow=2*square_size,ncol=square_size)
  splits<-floor(square_size/n)
  square_coords <- c()
  for (i in 1:n){
    if (i == n){
      coord <- square_size
    }
    else {
      coord <- i*splits
    }
    square_coords <- c(square_coords,coord)
    }
  pair_coords<-combn(square_coords,2)
  coord_pairs<-cbind(pair_coords,rbind(square_coords,square_coords))
  # Run the new subsets and collect them in the main matrix
  
  for (i in 1:dim(coord_pairs)[2]) { 
    if (coord_pairs[1,i] == coord_pairs[2,i]){
      temp_geno <- genotype[,(coord_pairs[1,i]-splits+1):coord_pairs[1,i]]
      temp_geno_rev <- 0
      if (model == 2){
        temp_geno_rev <- genotype_rev[,(coord_pairs[1,i]-splits+1):coord_pairs[1,i]]
      }
      coord_splits <-triangular_split_test(dim(temp_geno)[2],threads)
      partial_results <- foreach(j = 1:threads, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations_triangular(temp_geno,temp_geno_rev,phenotype=phenotype,coord_splits[j,],model=model) 
      return(subset)
      }
      snp_matrix[((coord_pairs[1,i]-splits+1)*2-1):(coord_pairs[1,i]*2),(coord_pairs[2,i]-splits+1):coord_pairs[2,i]] <- partial_results
    }
    else {
      if (square_size == coord_pairs[1,i]){
          temp_geno_1 <- genotype[,square_coords[n-1]:square_coords[n]]
          temp_geno_2 <- genotype[,(coord_pairs[2,i]-splits+1):coord_pairs[2,i]]
          temp_geno_rev <- 0
          if (model == 2)
          {
            temp_geno_rev <-  genotype_rev[,(coord_pairs[2,i]-splits+1):coord_pairs[2,i]]
          } 
          coords_splits <-square_split(threads,dim(temp_geno_1)[2])
          partial_results <- foreach(j = 1:threads, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
            subset <- partial_correlations(temp_geno_1[,coords_splits[j,1]:coords_splits[j,2]],temp_geno_2,temp_geno_rev,phenotype=phenotype,model=model)
            return(subset)
          }
         snp_matrix[(square_coords[n-1]*2-1):(square_coords[n]*2),(coord_pairs[2,i]-splits+1):coord_pairs[2,i]] <- partial_results
      }
      else if (square_size == coord_pairs[2,i]){ 
          temp_geno_2 <- genotype[,square_coords[n-1]:square_coords[n]]
          temp_geno_1 <- genotype[,(coord_pairs[1,i]-splits+1):coord_pairs[1,i]]
          print(dim(temp_geno_1))
          print(dim(temp_geno_2))
          temp_geno_rev <- 0
          if (model == 2){
            temp_geno_rev <- genotype_rev[,square_coords[n-1]:square_coords[n]]
          }
          coords_splits <-square_split(threads,dim(temp_geno_1)[2])
          partial_results <- foreach(j = 1:threads, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
            subset <- partial_correlations(temp_geno_1[,coords_splits[j,1]:coords_splits[j,2]],temp_geno_2,temp_geno_rev,phenotype=phenotype,model=model)
            return(subset)
          }
          #snp_matrix[((coord_pairs[1,i]-splits+1)*2-1):(coord_pairs[1,i]*2),square_coords[n-1]:square_coords[n]] <- partial_results
          }
      else { 
        temp_geno_1 <- genotype[,(coord_pairs[1,i]-splits+1):coord_pairs[1,i]]
        temp_geno_2 <- genotype[,(coord_pairs[2,i]-splits+1):coord_pairs[2,i]]
        print(dim(temp_geno_1))
        print(dim(temp_geno_2))
        temp_geno_rev <- 0
        if (model == 2){
          temp_geno_rev <- genotype_rev[,(coord_pairs[2,i]-splits+1):coord_pairs[2,i]]
        }
        coords_splits <-square_split(threads,dim(temp_geno_1)[2])
        partial_results <- foreach(j = 1:threads, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
          subset <- partial_correlations(temp_geno_1[,coords_splits[j,1]:coords_splits[j,2]],temp_geno_2,temp_geno_rev,phenotype=phenotype,model=model)
          return(subset)
        }
        #snp_matrix[((coord_pairs[1,i]-splits+1)*2-1):(coord_pairs[1,i]*2),(coord_pairs[2,i]-splits+1):coord_pairs[2,i]] <- partial_results
      }
    } 
  }
  }
  #return(snp_matrix)
}  

result<-memory_subspace(1,20,genotype = genotype,genotype_rev = genotype,phenotype=phenotype[,3],model = 1)

result<-memory_subspace(1,20,genotype = new_new_geno,genotype_rev = new_geno,phenotype=phenotype[,3],model = 1)

new_new_geno <- cbind(new_geno,new_geno)
dim(new_geno)

result<-epistatic.correlation_test(phenotype[,3],genotype,threads= 40,test=T,simple=T)

  # New partial correlation formula based on the non-traingular subsetting
partial_correlations <- function(genotype1,genotype2,genotype2_rev,phenotype,model=1){ 
    size1<-dim(genotype1)[2]
    size2<-dim(genotype2)[2]
    data_matrix <- matrix(0,nrow = 2*(size1),ncol=size2) 
    if (model==1){
      for (i in 1:size1){
          for (j in 1:size2){
            tmp_model = fastLm(phenotype ~ I(genotype1[,i])+I(genotype2[,j])+I(genotype1[,i]*genotype2[,j]))
            data_matrix[(i*2-1):(i*2),j]<-c(tmp_model$coefficients[length(tmp_model$coefficients)],summary(tmp_model)$coefficients[dim(summary(tmp_model)$coefficients)[1],4])
        }
      }
    }
    if (model==2){
      for (i in 1:size1){
        for (j in 1:size2){
          tmp_model = fastLm(phenotype ~ I(genotype1[,i])+I(genotype2[,j])+I(genotype1[,i]*genotype2_rev[,j]))
          data_matrix[(i*2-1):(i*2),j]<-c(tmp_model$coefficients[length(tmp_model$coefficients)],summary(tmp_model)$coefficients[dim(summary(tmp_model)$coefficients)[1],4])
        }
      }
    }
    gc(verbose = F)
    return(data_matrix)
  } 


partial_correlations_triangular <- function(genotype_1,genotype_rev_1,phenotype,coords,model=1){
  n=dim(genotype_1)[2]
  data_matrix <- matrix(0,nrow = 2*(coords[2]-coords[1]+1),ncol=dim(genotype_1)[2])
  matrix_row <- 0
  if (model==1){
    for (i in coords[1]:coords[2]){
      matrix_row <- matrix_row+1
      if (i < n){
        for (j in (i+1):n){
          tmp_model = fastLm(phenotype ~ I(genotype_1[,i])+I(genotype_1[,j])+I(genotype_1[,i]*genotype_1[,j]))
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
          tmp_model = fastLm(phenotype ~ I(genotype_1[,i])+I(genotype_1[,j])+I(genotype_1[,i]*genotype_rev_1[,j]))
          data_matrix[(matrix_row*2-1):(matrix_row*2),j]<-c(tmp_model$coefficients[length(tmp_model$coefficients)],summary(tmp_model)$coefficients[dim(summary(tmp_model)$coefficients)[1],4])
        }
      }
    }
  }
  return(data_matrix)
}

partial_correlations <- function(genotype1,genotype2,genotype2_rev,phenotype,model=1){ 
  size1<-dim(genotype1)[2]
  size2<-dim(genotype1)[2]
  genotype2 <- 0
  #print(c(size1,size2)) 
  data_matrix <- matrix(0,nrow = 2*(size1),ncol=size2) 
  if (model==1){
    for (i in 1:size1){
      for (j in 1:size2){
        tmp_model = fastLm(phenotype ~ I(genotype1[,i])+I(genotype1[,j])+I(genotype1[,i]*genotype1[,j]))
        data_matrix[(i*2-1):(i*2),j]<-c(tmp_model$coefficients[length(tmp_model$coefficients)],summary(tmp_model)$coefficients[dim(summary(tmp_model)$coefficients)[1],4])
      }
    }
  }
  if (model==2){
    for (i in 1:size1){
      for (j in 1:size2){
        tmp_model = fastLm(phenotype ~ I(genotype1[,i])+I(genotype2[,j])+I(genotype1[,i]*genotype2_rev[,j]))
        data_matrix[(i*2-1):(i*2),j]<-c(tmp_model$coefficients[length(tmp_model$coefficients)],summary(tmp_model)$coefficients[dim(summary(tmp_model)$coefficients)[1],4])
      }
    }
  }
  return(data_matrix)
} 


