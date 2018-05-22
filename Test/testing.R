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

genotype<-generate.genotype(ped = "/home/victor/Documents/ADD_files/INDICES_1_2.ped",tped = "/home/victor/Documents/ADD_files/INDICES1_1_2.tped",major.freq = 0.90)
genotype<-genotype[!(rownames(genotype) %in% ids_delete[,1]),]
tped <- fread("/home/victor/Documents/ADD_files/INDICES1_1_2.tped", data.table = F)




#### Developing Blocks ####

correlation_blocks <-function(genotype,threshold=0.9,measure="r2"){ 
  if (!(measure %in% c("r2","Dd"))){
    print("unknown LD measure, should be r2 or Dd")
  }
  else{
  snp_block_matrix <- as.matrix(rep(0,dim(genotype)[2]))
  rownames(snp_block_matrix) <- colnames(genotype)
  n_snp <- 1
  start <- 2
  n_block <- 1
  snps <- dim(genotype)[2]
  snp_block_matrix[n_snp,] <- n_block
  block_coords <- list()
  r_vector <- c()
  D_vector <- c()
  Dd_vector <- c()
  while(start <= snps){
      total <- sum(!is.na(genotype[,n_snp]+genotype[,start]))
      usable_values <- !is.na(genotype[,n_snp]+genotype[,start])
      p_AB <- (sum(genotype[usable_values,n_snp]+genotype[usable_values,start] == 0,na.rm = T)+(sum(genotype[usable_values,n_snp]+ genotype[usable_values,start] == 1 ,na.rm = T)/2)+(sum(genotype[usable_values,n_snp]*genotype[usable_values,start] == 1 ,na.rm = T)/2))/total
      p_A <- (sum(genotype[usable_values,n_snp] == 0,na.rm = T)+ sum(genotype[usable_values,n_snp] == 1,na.rm = T)/2)/total
      p_a <- 1-p_A
      p_B <- (sum(genotype[usable_values,start] == 0,na.rm = T)+ sum(genotype[usable_values,start] == 1,na.rm = T)/2)/total
      p_b <- 1-p_B
      D <- p_AB-(p_A*p_B)
      if( D < 0){
        D_min <- max(c(-(p_A*p_B),-(p_a*p_b)))
        Dd <- D/D_min
      }
      else {
        D_min <- min(c(p_A*p_b,p_a*p_B))
        Dd <- D/D_min
      }
      r2 <- D^2/(p_A*p_a*p_B*p_b)
      r_vector <- c(r_vector,r2)
      D_vector <- c(D_vector,D)
      Dd_vector <- c(Dd_vector,Dd)
      if (measure == "r2"){
        similarity <- r2
      }
      if (measure == "Dd"){
        similarity <- Dd
      }
      if (similarity >= threshold & similarity){ 
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
        }
      }
  }
  genotype <- as.matrix(genotype)
  new_genotype <- matrix(0L, nrow = dim(genotype)[1], ncol = n_block)
  counter <- 0
  tagging_markers <- c()
  message("Creating Blocks")
  for (coords in block_coords){
    counter <- counter + 1
    if (coords[1] == coords[2]){
      new_genotype[,counter] <- genotype[,coords[1]]
      tagging_markers <- c(tagging_markers,coords[1])
    }
    else {
        selection<- round(median(c(coords[1]:coords[2]))) 
        new_genotype[,counter] <- genotype[,selection] 
        tagging_markers <- c(tagging_markers,selection) 
      }  
    }
  output<-list(new_genotype,block_coords,snp_block_matrix,r_vector,D_vector,Dd_vector,tagging_markers)
  names(output)<-c("genotype","block_coords","snp_id_blocks","r_values","D_values","Dd_values","tagging_markers")
  return(output)
  }
}

genotype_cut<-correlation_blocks(genotype,threshold = 0.8,measure="r2")
geno_block<-(genotype_cut$genotype)
dim(geno_block)

geno_block[,2303]
blocks1<-genotype_cut[[1]]
blocks2<-genotype_cut[[2]]
blocks3<-genotype_cut[[3]]
blocks4<-genotype_cut[[4]]
blocks5<-genotype_cut[[5]]
blocks6<-genotype_cut[[6]]

max(blocks5)
min(blocks5)
hist(blocks4)
hist(blocks5)
hist(blocks6)

dim(genotype)

head(blocks2)

head(blocks3)


model



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
  
  
} 


blocks<-correlation_blocks_running(genotype)
tail(blocks)


correlation_blocks_network <-function(genotype,threshold=0.9,max_block_size=1000){ 
  snp_block_matrix <- as.matrix(rep(0,dim(genotype)[2]))
  rownames(snp_block_matrix) <- colnames(genotype)
  n_snp <- 1
  start <- 2
  n_block <- 1
  network_size <- 0
  r2_values <- c()
  block_coords <- list()
  snps <- dim(genotype)[2]
  snp_block_matrix[n_snp,] <- n_block
  while(start <= snps){ 
    if (n_snp == snps){
      snp_block_matrix[snps,1] = n_block
    }
    else{ 
      temp_sim <- 0
      network_pairs<-combn(c(n_snp:(start)),2)
      for (i in 1:dim(network_pairs)[2]){
        start_t <- network_pairs[1,i]
        n_snp_t <- network_pairs[2,i]
        # if (n_snp == snps){
        #   snp_block_matrix[snps,1] = n_block
        # }
        #else{
        #matches<-sum(genotype[,n_snp]/genotype[,start] == 1,na.rm = T)+(sum(c(c(genotype[,n_snp] == 1.5)+ c(genotype[,start] == 1.5)) == 1 ,na.rm = T)/2)
        total <- sum(!is.na(genotype[,n_snp_t]+genotype[,start_t]))
        usable_values <- !is.na(genotype[,n_snp_t]+genotype[,start_t])
        #similarity <- matches/total
        p_AB <- (sum(genotype[usable_values,n_snp_t]+genotype[usable_values,start_t] == 0,na.rm = T)+(sum(genotype[usable_values,n_snp_t]+ genotype[usable_values,start_t] == 1 ,na.rm = T)/2)+(sum(genotype[usable_values,n_snp_t]*genotype[usable_values,start_t] == 1 ,na.rm = T)/2))/total
        p_A <- (sum(genotype[usable_values,n_snp_t] == 0,na.rm = T)+ sum(genotype[usable_values,n_snp_t] == 1,na.rm = T)/2)/total
        p_a <- 1-p_A
        p_B <- (sum(genotype[usable_values,start_t] == 0,na.rm = T)+ sum(genotype[usable_values,start_t] == 1,na.rm = T)/2)/total
        p_b <- 1-p_B
        #print(c(p_A,p_B,p_AB,start))
        D <- p_AB-(p_A*p_B)
        r2 <- D^2/(p_A*p_a*p_B*p_b)
        r2_values <- c(r2_values,r2)
        temp_sim <- r2+temp_sim
        network_size <- network_size + 1
      }
      similarity <- temp_sim
      print(c(similarity,similarity/network_size,network_size,start,n_snp))
      if (similarity/network_size >= threshold && network_size < 1001){
        snp_block_matrix[start,1] = n_block
        start <- start+1
        network_size <-0
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
        n_block <- n_block +1
        snp_block_matrix[start,1] = n_block
        start <- start +1
        network_size <- 0
        }
      }
    }
  }
  genotype <- as.matrix(genotype)
  new_genotype <- matrix(0L, nrow = dim(genotype)[1], ncol = n_block)
  counter <- 0
  tagging_markers <- c()
  message("Creating Blocks")
  for (coords in block_coords){
    counter <- counter + 1
    if (coords[1] == coords[2]){
      new_genotype[,counter] <- genotype[,coords[1]]
      tagging_markers <- c(tagging_markers,coords[1])
    }
    else {
      selection<- round(median(c(coords[1]:coords[2]))) 
      new_genotype[,counter] <- genotype[,selection] 
      tagging_markers <- c(tagging_markers,selection) 
    }  
  }
  output<-list(new_genotype,tagging_markers,snp_block_matrix,r2_values)
  names(output)<-c("genotype","tagging_markers","genotype_block_matrix")
  return(output)
}

blocks<-correlation_blocks_network(genotype,0.8)


max(table(blocks$genotype_block_matrix))

hist(blocks$r2)

dim(blocks[[1]])
dim(blocks[[2]])
tail(blocks[[2]])
tail(blocks[[1]])

length(blocks$tagging_markers)
blocks<-blocks$genotype_block_matrix
max(table(blocks))


max(hist(blocks[[2]]))
dim(blocks[[1]])
tail(blocks[[3]])



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
genotype<-generate.genotype(ped = "/data/nbg153/ADHD/WISH_files_newest/trimmed_ped",tped = "/data/nbg153/ADHD/WISH_files_newest/trimmed_tped",major.freq = 0.9,coding="linear")
genotype<-genotype[!(rownames(genotype) %in% ids_delete[,1]),]
tped <- fread("/data/nbg153/ADHD/WISH_files_newest/INDICES1_1_2.tped", data.table = F)

dim(genotype)
ped <- fread("/data/nbg153/ADHD/WISH_files_newest/trimmed_ped", data.table = F)
tped <- fread("/data/nbg153/ADHD/WISH_files_newest/trimmed_tped", data.table = F)



genotype<-generate.genotype(ped = ped,tped = "/data/nbg153/ADHD/WISH_files_newest/trimmed_tped",major.freq = 0.9)

new_geno <- cbind(genotype,genotype,genotype)

dim(genotype)

system.time(result<-epistatic.correlation(phenotype[,3],genotype,threads= 40,test=F,simple=F))
result2<-epistatic.correlation(phenotype[,3],geno_block,threads= 40,test=F,simple=F)

hist(result2$Pvalues)
r2_bf<-p.adjust(result2$Pvalues)
r2_fdr<-p.adjust(result2$Pvalues,method="fdr")
min(r2_bf)
min(r2_fdr)

min(result2$Pvalues)
which(result$Pvalues < 5.950226e-08)


3204039%%2695
names<-row.names(result$Pvalues)

tped[tped$V2==names[2379],1:5]
tped[tped$V2==names[1189],1:5]

which.min(result$Pvalues[2379,1189])

result$Pvalues[,2379]

dim(result$Pvalues)
r2_bf<-p.adjust(result$Pvalues)
r2_fdr<-p.adjust(result$Pvalues,method="fdr")
min(r2_bf)
min(r2_fdr)
which(r2_fdr < 0.05)

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
            tmp_model = glm(phenotype ~ I(genotype1[,i])+I(genotype2[,j])+I(genotype1[,i]*genotype2[,j]))
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



genotype_new<- genotype
genotype_new[genotype_new==2] <- 0
genotype_new[genotype_new==1] <- 4
genotype_new[genotype_new==1.5] <- 1
genotype_new[genotype_new==4] <- 2

values <- c()
for (i in 1:100){
  values <- c(values,max(rnorm(i)))
}
  
plot(1:100,values)
  
  
qnorm(0.1) 
?qnorm




# GLM testing

partial_correlations_triangular_glm <- function(genotype_1,genotype_rev_1,phenotype,coords,model=1){
  n=dim(genotype_1)[2]
  data_matrix <- matrix(0,nrow = 2*(coords[2]-coords[1]+1),ncol=dim(genotype_1)[2])
  matrix_row <- 0
  if (model==1){
    for (i in coords[1]:coords[2]){
      matrix_row <- matrix_row+1
      if (i < n){
        for (j in (i+1):n){
          tmp_model = glm(phenotype ~ I(genotype_1[,i])+I(genotype_1[,j])+I(genotype_1[,i]*genotype_1[,j]),family=binomial())
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
          tmp_model = glm(phenotype ~ I(genotype_1[,i])+I(genotype_1[,j])+I(genotype_1[,i]*genotype_rev_1[,j]),family=binomial())
          data_matrix[(matrix_row*2-1):(matrix_row*2),j]<-c(tmp_model$coefficients[length(tmp_model$coefficients)],summary(tmp_model)$coefficients[dim(summary(tmp_model)$coefficients)[1],4])
        }
      }
    }
  }
  return(data_matrix)
} 


partial_correlations_glm <- function(genotype_1,genotype_rev_1,phenotype,coords,model=1){
  n=dim(genotype_1)[2]
  data_matrix <- matrix(0,nrow = 2*(coords[2]-coords[1]+1),ncol=dim(genotype_1)[2])
  matrix_row <- 0
  if (model==1){
    for (i in coords[1]:coords[2]){
      matrix_row <- matrix_row+1
      if (i < n){
        for (j in (i+1):n){
          tmp_model = glm(phenotype ~ I(genotype_1[,i])+I(genotype_1[,j])+I(genotype_1[,i]*genotype_1[,j]),family=binomial())
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
          tmp_model = glm(phenotype ~ I(genotype_1[,i])+I(genotype_1[,j])+I(genotype_1[,i]*genotype_rev_1[,j]),family=binomial())
          data_matrix[(matrix_row*2-1):(matrix_row*2),j]<-c(tmp_model$coefficients[length(tmp_model$coefficients)],summary(tmp_model)$coefficients[dim(summary(tmp_model)$coefficients)[1],4])
        }
      }
    }
  }
  return(data_matrix)
} 

epistatic.correlation_glm <- function(phenotype,genotype,threads=1,test=T,simple=T){ 
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
    decide_1<-(genotype_rev==0)
    decide_2<-(genotype_rev==2)
    genotype_rev[decide_1] <- 2
    genotype_rev[decide_2] <- 0
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
      subset <- partial_correlations_glm(genotype[,1:316],genotype_rev[,1:316],phenotype=phenotype,coord_splits[j,],model=model) 
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
      subset <- partial_correlations_glm(genotype[,1:316],genotype_rev[,1:316],phenotype=phenotype,coord_splits[j,],model=2) 
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
    coord_splits <-glm_split(dim(genotype)[2],threads)
    snp_matrix <- foreach(j = 1:threads, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations_glm(genotype,genotype_rev,phenotype=phenotype,coord_splits[j,],model=model) 
      return(subset)
    }
    # Running opposite minor/major co-linearity model
    coord_splits <-triangular_split(dim(genotype)[2],threads)
    snp_matrix_rev <- foreach(j = 1:threads, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations_glm(genotype,genotype_rev,phenotype=phenotype,coord_splits[j,],model=model) 
      return(subset)
    }
  }
  else if (test==F && simple==F) {
    coord_splits <-triangular_split(dim(genotype)[2],threads)
    snp_matrix <- foreach(j = 1:threads, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations_glm(genotype,genotype_rev,phenotype=phenotype,coord_splits[j,],model=model) 
      return(subset)
    }
    # Running opposite minor/major co-linearity model
    coord_splits <-triangular_split(dim(genotype)[2],threads)
    snp_matrix_rev <- foreach(j = 1:threads, .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      subset <- partial_correlations_glm(genotype,genotype_rev,phenotype=phenotype,coord_splits[j,],model=model) 
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
      subset <- partial_correlations_glm(genotype,genotype_rev,phenotype=phenotype,coord_splits[j,],model=model) 
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

rand_pheno <- sample(c(1,0),193,T)
system.time(result2<-epistatic.correlation(rand_pheno,genotype[,1:1000],threads= 10,test=T,glm=T,simple=F)) 

tmp_model<-speedglm(rand_pheno ~ I(genotype[,1])+I(genotype[,2])+I(genotype[,2]*genotype[,1]),family=binomial())


mode1<-formula(phenotype[,3] ~ I(genotype[,1])+I(genotype[,2])+I(genotype[,2]*genotype[,1]))
mode2<-formula(phenotype[,3] ~ I(genotype[,1])+I(genotype[,2])+I(genotype[,2]*genotype[,1]))

#visualization testing



effect_sumarize<-colSums(abs(e_e))
effect_sumarize<-abs(colSums(e_e))

plot(1:7758,effect_sumarize)



mins_p<-apply(p_e,2,min)
plot(1:7758,-log(mins_p))
hist(mins_p)

which.min(mins_p)

effect_sumarize<-colSums(-log(p_e))
hist(effect_sumarize)
which(effect_sumarize > 18000)
hist(effect_sumarize)
effect_sumarize<-abs(colSums(p_e))

plot(1:7758,effect_sumarize)

tped[tped[,2] == "BovineHD1300004594",]

table(tped[,1])
effect_test<-c()
for (i in 1:8000){
  print(i)
  effect_test <- c(effect_test,sum(-log(runif(1000))))
}

hist(effect_test) 






 
#' Function for summarizing individual variant epistasis results 
#' @description Extract summary information from the epistasis analysis 
#' for each variant
#' @usage variant_summary(tped,correlations)
#' @param tped Input tped file as .tped file or data frame. The tped file (.tped)
#' is a transposed ped file, from Plink. 
#' This file contains the SNP and genotype information where one row is a SNP.
#' The first 4 columns of a TPED file are the same as a 4-column MAP file.
#' Then all genotypes are listed for all individuals for each particular SNP on 
#' each line. Again, SNPs are 1,2-coded.#'
#' @param correlations List of epistatic correlations and p-values generated by
#' epistatic.correlation()
#' @return Plots a pseudo manhattan plot
#' @examples
#' pseudo_manahattan(tped,correlations)
#' 
#' @export


summary_matrix <- function(tped,correlations) {
  matrix_base <- tped[tped[,2] %in% rownames(correlations$Pvalues),c(1,2,4)]
  likelyhood_sum <- colSums(-log(correlations$Pvalues))
  effect_sum <- colSums(abs(correlations$Coefficients))
  min_p <- apply(correlations$Pvalues,2,min)
  summary_matrix<-cbind(matrix_base,likelyhood_sum,effect_sum,min_p)
  colnames(summary_matrix) <- c("chr","variant","position","likelihood_sum","Effect_Sum","min_Pvalue")
  return(summary_matrix)
} 



pseudo_manhattan(tped,model_GIFT,values="c") 
test<-variant_matrix(tped,model_GIFT)

tped<-data.table::fread("../caw_project/GIFT_final2.tped",data.table = F)

tped <- tped[!(tped$V1 %in%  c(0,23,24,25,26,27,28,29,30,31,32,33)),]

load("../../nbg153/caw_project/model_GIFT.RData")
genome.interaction(tped,model_GIFT)

correlations <- model_GIFT
genome.interaction <- function(tped,correlations,quantile=0.9) {
  if (is.character(tped)){
    message("loading tped file")
    tped <- fread(tped,data.table=F)
  }
  else if (!is.data.frame(tped)){
    stop("tped file not file or data frame")
  }
  new_P <- (1-correlations$Pvalues)
  map<-tped[tped[,2] %in% rownames(correlations$Pvalues),1:2]
  counter = 0
  ends <- c()
  chr_list <- c()
  for (i in map[,1]){
    if (counter == 0){
      counter <- counter + 1
      starts <- 1
      chr <- i
      chr_list <- c(chr_list,chr)
    }
    else {
      if (chr == i){
        counter <- counter + 1 
        if (counter == dim(map)[1]){
          ends <- c(ends,counter)
        }
      }
      if (chr != i){
        chr <- i
        chr_list <- c(chr_list,chr)
        ends <- c(ends,counter)
        counter <- counter + 1 
        starts <- c(starts,counter)
      }
    }
  }
  coord_splits<-cbind(starts,ends)
  visualization_matrix <- matrix(nrow = length(starts),ncol = length(starts))
  colnames(visualization_matrix) <- chr_list
  rownames(visualization_matrix) <- chr_list
  for (i in 1:length(starts)){
    for (j in 1:length(starts)){
      subset <- c(new_P[coord_splits[i,1]:coord_splits[i,2],coord_splits[j,1]:coord_splits[j,2]])
      subset <- abs(subset)
      visualization_matrix[i,j] <- quantile(subset,quantile,na.rm=T)
    }
  }
  visualization_matrix <- 2*(visualization_matrix-min(visualization_matrix))/(max(visualization_matrix)-min(visualization_matrix))-1
  corrplot(visualization_matrix, type="upper",title= "Pairwise Chromosome Interaction Map",mar=c(0,0,2,0))
}


genome.interaction(tped,model_GIFT)

hist(test$likelihood_sum)
hist(test$Effect_Sum)
hist(test$chr)

dim(correlations$Pvalues)


sim_genotype <- function(samples,snps) {
  sp<-sample(c(0,1,2),samples*snps,T)  
  genotype_matrix <- matrix(sp,nrow=samples,ncol=snps)
  return(genotype_matrix)
}
i=3000
for (i in c(100,500,1000,2000,4000)){
  genotype<-sim_genotype(i,2000)
  phenotypes <- rnorm(i)
  time<-system.time(epistatic.correlation(phenotypes,genotype,40,simple = T,test=F))
  print(time)
}


for (i in c(5,10,15,20,30,40)){
  genotype<-sim_genotype(500,2000)
  phenotypes <- rnorm(500)
  time<-system.time(epistatic.correlation(phenotypes,genotype,simple = T,test=F))
  print(time)
}

for (i in c(5,10,15,20,30,40)){
  genotype<-sim_genotype(500,3000)
  phenotypes <- rnorm(500)
  time<-system.time(epistatic.correlation(phenotypes,threads = 40,genotype,simple = T,test=F))
  print(time)
}



run_1000 <- c(261,140,104,86,58,53)
run_2000 <- c(1038,559,395,327,231,206)
run_3000 <- c(4543,2463,1757,1441,1113,953)/2
cores <- c(5,10,15,20,30,40)




run_data<-as.data.frame(cbind(c(cores,cores,cores),c(run_1000,run_2000,run_3000),c(1000,1000,1000,1000,1000,1000,2000,2000,2000,2000,2000,2000,3000,3000,3000,3000,3000,3000)))

colnames(run_data) <- c("cores","runtime","size")

run_data$size <- as.character(run_data$size)

library(ggplot2)

gg <- ggplot(run_data,aes(x=cores, y=(runtime),col=size )) + 
  geom_point()+
  ggtitle("Multi-thread Scaling")+
  ylab("Seconds")+
  xlab("Threads")
gg <- gg+guides(size=F)
gg <- gg+labs(col="N-variants")
gg

seconds <- c(174,216,277,346,466,551)
samples <- c(100,500,1000,2000,3000,4000)
sample_runtime<-as.data.frame(cbind(seconds,samples))

gg<-ggplot(sample_runtime,aes(x=samples, y=seconds)) + 
  geom_point(colour="blue",size=4)+
  ggtitle("Sample Size Scaling")+
  ylab("seconds")+
  xlab("N-samples")
gg<-gg+guides(size=FALSE)
gg


ygenotype<-sim_genotype(500,1000)
phenotypes <- rnorm(500)
epistatic.correlation(phenotypes,threads = 40,genotype,simple = F,test=F)

class(tped)

genome.interaction(tped,model_GIFT)

pairwise.chr.map <- function(chr1,chr2,tped,correlations,span=10^6) {  
  new_P <- (1-correlations$Pvalues)
  message("loading tped file")
  tped <- fread(tped,data.table=F)
  total_map <- tped[tped[,2] %in% rownames(correlations$Coefficients),c(1,4)]
  total_map[,2] <- as.numeric(total_map[,2])
  snps1<-c(which(total_map[,1] == chr1))
  snps2 <-c(which(total_map[,1] == chr2))
  values <- abs(correlations$Pvalues[snps1,snps2])
  #values <- -log(values)
  #threshold<-quantile(values,0.95)
  #values[values < threshold] <- 0
  
  heatmap3(values,Rowv = NA,Colv = NA)
}
par(oma=c(0,0,2,0))
pairwise.chr.map(1,2,"../caw_project/GIFT_final2.tped",model_GIFT,grid_size = 50)


hist(model_GIFT$Coefficients)

#heatmap3(visualization_matrix,scale="none",main="Pairwise Chromosomal Interaction",Rowv = NA,Colv = NA,xlab=xlabel,ylab=ylabel ,labRow=c("start",rep("",dim(chromosome_choords1)[1]-2),"end"),labCol=c("start",rep("",dim(chromosome_choords2)[1]-2),"end"))







#' Function to plot summary pseudo-manhattan plots of variants.
#' @description Visualize summary statistics for interactions based on total sum 
#' of -loglikelihoods for each variant across all interactions og the sum of
#' effect sizes.
#' @import ggplot2 
#' @usage pseudo_manahattan(tped,correlations)
#' @param tped Input tped file as .tped file or data frame. The tped file (.tped)
#' is a transposed ped file, from Plink. 
#' This file contains the SNP and genotype information where one row is a SNP.
#' The first 4 columns of a TPED file are the same as a 4-column MAP file.
#' Then all genotypes are listed for all individuals for each particular SNP on 
#' each line. Again, SNPs are 1,2-coded.#'
#' @param correlations List of epistatic correlations and p-values generated by
#' epistatic.correlation()
#' @return Plots a pseudo manhattan plot
#' @examples
#' pseudo_manahattan(tped,correlations)
#' 
#' @export


pseudo_manhattan<- function(tped,correlations,values="p"){
  if (is.character(tped)){
    message("loading tped file")
    tped <- fread(tped,data.table=F)
  }
  map<-tped[tped[,2] %in% rownames(correlations$Pvalues),1:2]
  if ((sum(map[,2] == rownames(correlations$Pvalues)))== dim(correlations$Pvalues)[1]){
    if (values=="p"){
      likelyhood_sum <- colSums(-log(correlations$Pvalues))
      #plot(1:length(likelyhood_sum),likelyhood_sum,col=map[,1]%%2+3,xlab="N-variant",ylab="Sum of log-likelihood",main="Pseudo-Manhattan Plot")
      data_p<-as.data.frame(cbind(likelyhood_sum,c(1:length(likelyhood_sum)),map[,1]))
      colnames(data_p)<-c("like_sum","Nvar","chr")
      data_p$chr <- as.factor(data_p$chr)
      gg<-ggplot(data_p,aes(x=Nvar, y=like_sum,col=chr ))+ geom_point()+
        ggtitle("Pseudo-Manhattan Plot")+
        xlab("N-variant")+
        ylab("Sum of -log-likelihood")
      gg <- gg+guides(col=F)
      return(gg)
    }
    if (values=="c"){
      likelyhood_sum <- colSums(abs(correlations$Coefficients))
      data_p<-as.data.frame(cbind(likelyhood_sum,c(1:length(likelyhood_sum)),map[,1]))
      colnames(data_p)<-c("like_sum","Nvar","chr")
      data_p$chr <- as.factor(data_p$chr)
      gg<-ggplot(data_p,aes(x=Nvar, y=like_sum,col=chr ))+ geom_point()+
        ggtitle("Pseudo-Manhattan Plot")+
        xlab("N-variant")+
        ylab("Sum of Effect Sizes")
      gg <- gg+guides(col=F)
      return(gg)
    }
  }
  else{
    message("tped and model output variant names do not match")
  }
}

#' * genotype The tagging genotypes selected from the blocks
#' * tagging_genotype The genotype selected to represent each block. The median genotype, rounded down is selected
#' * genotype_block_matrix A matrix indicating which block each genotype belongs to


pig_fe <-read.table("/data/nbg153/FeedOmics/phenotypes.txt")




ped <- fread("/data/nbg153/ADHD/WISH_files_newest/trimmed_ped", data.table = F)
tped <- fread("/data/nbg153/ADHD/WISH_files_newest/trimmed_tped", data.table = F)

dim(ped)

ped_small <- ped[,1:5006] 
tped_small <- tped[1:2500,]


names <- c()
for (i in 1:dim(ped)[1]){
  names <- c(names, paste("sample",as.character(i),sep=""))
}
names


snp_names <- c()
for (i in 1:dim(tped_small)[1]){
  snp_names <- c(snp_names, paste("snp",as.character(i),sep=""))
}
snp_names



tped_small[,2] <- snp_names
ped_small[,2] <- names


phenotype <-round(runif(203,1,100)) 

data<-as.data.frame(cbind(names,phenotype))

write.table(ped_small,"/data/nbg153/ADHD/WISH_files_newest/test.ped",col.names=F,row.names=F,quote=F)
write.table(tped_small,"/data/nbg153/ADHD/WISH_files_newest/test.tped",col.names=F,row.names=F,quote=F)
write.table(data,"/data/nbg153/ADHD/WISH_files_newest/test_pheno.txt",col.names=F,row.names=F,quote=F)

chr1 <- 1
chr2 <- 2
region1 <- 1


correlations<-model_GIFT  


tped[,1]<- chr

genotype<-generate.genotype(ped,tped)



LD_genotype<-LD_blocks_t(genotype)
LD_genotype$tagging_genotype
genotype1 <- LD_genotype$genotype

results<-epistatic.correlation(phenotype[,2], genotype,threads = 20 ,test=F)

phenotype[,2] <- rnorm(203)+sample(c(1:10),203,replace=T)

write.table(phenotype,"test_pheno.txt",col.names=F,row.names=F,quote=F)

correlations<-epistatic.correlation(phenotype[,2], genotype, threads = 20 ,test=F,simple=F)
hist(correlations$Coefficients)


genome.interaction(tped,correlations) 

write.table(tped,"tped.test",col.names=F,row.names=F,quote=F)
correlations$Coefficients[correlations$Coefficients > 1000 |  correlations$Coefficients < -1000] <-0
hist(correlations$Coefficients)

modules <-generate.modules(correlations,threads=2)


generate.modules <- function(correlations,values="Coefficients",power=c(seq(1,10,0.1),c(12:22)),n.snps=dim(correlations$Coefficients)[1],minClusterSize=50,type="unsigned",threads=2) {
  enableWGCNAThreads(threads)
  if (values=="Pvalue"){
    corr <- (1-correlations$Pvalues)*(correlations$Coefficients/abs(correlations$Coefficients))
  }
  if (values=="Coefficients"){
    temp_corr<-correlations$Coefficients
    temp_corr[temp_corr < 0] <- 0
    temp_corr <- temp_corr/(max(temp_corr))
    corr <- temp_corr
    temp_corr <- correlations$Coefficients
    temp_corr[temp_corr > 0] <- 0
    temp_corr <- temp_corr/(abs(min(temp_corr)))
    corr[temp_corr < 0] <- temp_corr[temp_corr < 0]
  }
  sft = pickSoftThreshold(corr, powerVector = c(seq(1,10,0.1),c(12:30)), verbose = 5)
  connectivity <- adjacency.fromSimilarity(corr, power=sft$powerEstimate,type=type)
  sizeGrWindow(10,5)
  par(mfrow=c(1,2))
  hist(connectivity,xlab="connectivity")
  scaleFreePlot(connectivity)
  par(mfrow=c(1,1))
  #select SNPs for network construction based on connectivity
  select.snps <- corr[rank(-colSums(connectivity),ties.method="first")<=n.snps,rank(-colSums(connectivity),ties.method="first")<=n.snps ]
  select.snps[,c(1:ncol(select.snps))] <- sapply(select.snps[,c(1:ncol(select.snps))], as.numeric)
  #create adjacency matrix (correlation matrix raised to power beta)
  adjMat <- adjacency.fromSimilarity(select.snps,power=sft$powerEstimate,type=type)
  #calculate dissimilarity TOM
  dissTOM <- 1-(TOMsimilarity(adjMat))
  #create gene dendrogram based on diss TOM
  genetree <- flashClust(as.dist(dissTOM), method="average")
  #cut branches of the tree= modules
  dynamicMods = cutreeDynamic(dendro=genetree, distM=dissTOM, 
                              deepSplit=2,pamRespectsDendro=F, minClusterSize=minClusterSize)
  #give modules a color as name
  moduleColors = labels2colors(dynamicMods)
  #sizeGrWindow(8,6)
  #plot dendrogram with the module colors
  plotDendroAndColors(genetree, moduleColors,
                      dendroLabels=F, hang=0.03,
                      addGuide=T, guideHang = 0.05,
                      main="Gene dendrogram and modules")
  output <- list(select.snps,connectivity,adjMat,dissTOM,genetree,dynamicMods,moduleColors,sft$powerEstimate)
  names(output) <- c("SNPs","connectivity","adjMat","dissTom","genetree","modules","modulecolors","power.estimate")
  return(output)
}




#' Visualization of chromosomal pairwise region epistatic interaction strength, based on 
#' statistical significance 
#' @description Visualization of chromosome pairwise region epistatic interaction strength, based on 
#' statistical significance. The value is based of the most signficant epistatic interaction in each
#' region pair, ranging from 1 ( strongest) to 0 (weakest). By defaulty chromosomes are separated into
#' 1 Mb regions, but if SNPs are more spaced out that this it will adjust to the smallest region that fit
#' the data.  
#' @import heatmap3
#' @usage pairwise.chr.map(chr1,chr2,tped,correlations,span=10^6)
#' @param chr1 The name of the first chromosome in the comparison, matching the name
#' from the tped file
#' @param chr2 The name of the second chromosome in the comparison, matching the name
#' from the tped file
#' @param tped The tped file used in generate.genotype(). The SNPs must
#' be sorted by chromosome and position on the chromosome, matching the order of the SNPs in the correlation 
#' matrices. 
#' @param span Region in bp. Default is 1 Mb (10^6)
#' @param correlations List of epistatic correlations and p-values genrated by
#' epistatic.correlation()
#' @return Outputs a plot visualizing the pairwise chromosome region interaction
#' @examples
#'  pairwise.chr.map("1","2",tped,correlations)
#' 
#' @export


test_pairwise.chr.map <- function(chr1,chr2,tped,correlations,span=10^6) {  
  new_P <- (1-correlations$Pvalues)
  #message("loading tped file")
  #tped <- fread(tped,data.table=F)
  total_map <- tped[tped[,2] %in% rownames(correlations$Pvalues),c(1,4)]
  total_map[,2] <- as.numeric(total_map[,2])
  map <-total_map[total_map[,1] == chr1,]
  first_snp <- map[1,2]
  last_snp <- map[dim(map)[1],2]
  #size <- round((last_snp-first_snp)/span,digits=0)
  progress <- 1
  ends <- c()
  for (snp in map[,2]){   
    if ( snp < (first_snp+progress*span)) {
      starts <- 1
      row <- 1
    }
    else {
      while (snp > first_snp+(progress+1)*span) {
        progress <- progress +1
      }
      print("what")
      ends <- c(ends,row)
      row <- row +1
      starts <- c(starts,row)
      progress <- progress + 1
    }
  }
  ends<- c(ends,dim(map)[1])
  chromosome_choords1 <- cbind(starts,ends)
  chromosome_choords1 <- chromosome_choords1 + which(total_map[,1] == chr1)[1]-1
  map <-total_map[total_map[,1] == chr2,]
  first_snp <- map[1,2]
  last_snp <- map[dim(map)[1],2]
  progress <- 1
  ends <- c()
  for (snp in map[,2]){   
    if ( snp < (first_snp+progress*10^6)) {
      #print("!t")
      starts <- 1
      row <- 1
    }
    else {
      while (snp > first_snp+(progress+1)*10^6) {
        progress <- progress +1
      }
      #print("what")
      ends <- c(ends,row)
      row <- row +1
      starts <- c(starts,row)
      progress <- progress + 1
    }
  }
  ends<- c(ends,dim(map)[1])
  chromosome_choords2 <- cbind(starts,ends)
  chromosome_choords2 <- chromosome_choords2 + which(total_map[,1] == chr2)[1]-1
  visualization_matrix <- matrix(nrow = dim(chromosome_choords1)[1],ncol = dim(chromosome_choords2)[1])
  colnames(visualization_matrix) <- 1:(dim(chromosome_choords2)[1])
  rownames(visualization_matrix) <- 1:(dim(chromosome_choords1)[1]) 
  for (i in 1:dim(chromosome_choords1)[1]){
    for (j in 1:dim(chromosome_choords2)[1]){
      subset <- c(new_P[chromosome_choords1[i,1]:chromosome_choords1[i,2],chromosome_choords2[j,1]:chromosome_choords2[j,2]])
      subset <- abs(subset)
      visualization_matrix[i,j] <- max(subset)
      
    }
  }
  xlabel<-paste("Chromosome=",as.character(chr2),", N-regions=",as.character(dim(chromosome_choords2)[1]))
  ylabel<-paste("Chromosome=",as.character(chr1),", N-regions=",as.character(dim(chromosome_choords1)[1]))
  heatmap3(visualization_matrix,scale="none",main="Pairwise Chromosomal Interaction",Rowv = NA,Colv = NA,xlab=xlabel,ylab=ylabel ,labRow=c("start",rep("",dim(chromosome_choords1)[1]-2),"end"),labCol=c("start",rep("",dim(chromosome_choords2)[1]-2),"end"))
}

load("../caw_project/model_GIF05.RData")

tped <- fread("../caw_project/GIFT_final2_filtered05.tped",data.table = F)

p_e <- model_GIFT$Pvalues 
e_e <- model_GIFT$Coefficients 


correlations <- model_GIFT

pvalues <- correlations$Pvalues

chr1 <- 1
chr2 <- 2
grid_size <- 50


test_pairwise.chr.map <- function(chr1,chr2,tped,correlations,grid_size){
  total_map <- tped[tped[,2] %in% rownames(correlations$Pvalues),c(1,4)]
  chr_pair_values <-correlations$Pvalues[which(total_map[,1] == chr1),which(total_map[,1] == chr2)]
  chr_sizes <-dim(chr_pair_values)   
  chr1_coords<-1:chr_sizes[1]
  chr2_coords<-1:chr_sizes[2]
  chr1_chunk<-(length(chr1_coords)/grid_size)
  chr2_chunk<-floor(length(chr2_coords)/grid_size)
  chr1_splits<-split(chr1_coords,ceiling(seq_along(chr1_coords)/chr1_chunk))
  chr2_splits<-split(chr2_coords,ceiling(seq_along(chr2_coords)/chr2_chunk))
  mean_matrix<-matrix(0,nrow=length(chr1_splits),ncol=length(chr2_splits))
  counter <- c(0,0)
  for ( i in chr1_splits){
    counter[1]<- counter[1]+1
    counter[2]<-0
    for (j in chr2_splits){
      counter[2]<- counter[2]+1
           mean_matrix[counter[1],counter[2]] <- mean(chr_pair_values[i,j])
    }
  }
  xlabel<-paste("Chromosome=",as.character(chr2),", N-regions=",as.character(counter[1]))
  ylabel<-paste("Chromosome=",as.character(chr1),", N-regions=",as.character(counter[2]))
  heatmap3(mean_matrix,scale="none",main="Pairwise Chromosomal Interaction",Rowv = NA,Colv = NA,xlab=xlabel,ylab=ylabel,labRow=c("start",rep("",(counter[1]-2)),"end"),labCol=c("start",rep("",(counter[2]-2)),"end"))
}

test_pairwise.chr.map(1,2,tped,correlations,40) 



#### Simulate Genotype ####


uniform_genotype_sampler <- function(n_samples,n_genotypes) {
  genotypes <- matrix(0,nrow = n_samples, ncol= n_genotypes)
  for (i in c(1:n_genotypes)){
    maf<-runif(1,0,0.5)
    maf2 <- maf^2
    genotype <- runif(n_samples)
    mm <-genotype < maf2
    Mm <- genotype >= maf2 & genotype < (maf2+maf*(1-maf))
    MM <- genotype  >=(maf2+maf*(1-maf))
    genotype[mm] <- 2
    genotype[Mm] <- 1.5
    genotype[MM] <- 1
    genotypes[,i] <- genotype
  }
  return(genotypes)
}


blockwise_genotype_sampler <- function(n_samples,n_genotypes) {
  genotypes <- matrix(0,nrow = n_samples, ncol= n_genotypes)
  for (i in c(1:n_genotypes)){
    maf<-runif(1,0,0.5)
    maf2 <- maf^2
    genotype <- runif(n_samples)
    mm <-genotype < maf2
    Mm <- genotype >= maf2 & genotype < (maf2+maf*(1-maf))
    MM <- genotype  >=(maf2+maf*(1-maf))
    genotype[mm] <- 2
    genotype[Mm] <- 1.5
    genotype[MM] <- 1
    genotypes[,i] <- genotype
  }
  return(genotypes)
}





epistasis_layer <- function(genotype,n_pheno){
  n_snps <- dim(genotype)[2]
  epi_pheno <- rep(0,n_pheno)
  effect_sizes<- matrix(0,nrow=n_snps,ncol=n_snps) 
  for ( i in 1:n_snps){
    n_values <- length(i:n_snps)
    x1 <- 1000
    x0 <- 0.001
    y <- runif(n_values)
    n=5
    #power law
    x = ((x1^(n+1) - x0^(n+1))*y + x0^(n+1))^(1/(n+1))
    #rescaling to frequent low values
  
    effect_sizes[i,i:n_snps] <- x
  }
  x_max = max(effect_sizes)
  x <- effect_sizes
  x = (-1*x+max(x))
  x <- x^(log10(dim(genotype)[2]^2))
  x=x/max(x)*10
  for (i in 1:n_snps){
    sign <-sample(c(-1,1),replace=T,n_snps)
    x[,i] <- sign*x[,i]
  }
  x = x*(effect_sizes !=0)
  effect_sizes = x
  for (i in 1:n_snps){
    if (i+1 <= n_snps){
    for (j in (i+1):n_snps){
      epi_geno<-genotype[,i]*genotype[,j]      
      epi_pheno <- epi_pheno+effect_sizes[i,j]*epi_geno
    }  
    }
  }
 return(list(epi_pheno,effect_sizes)) 
}


main_layer <- function(genotype,range){ 
  x1 <- 1000
  x0 <- 0.001
  y <- runif(dim(genotype)[2])
  n=5
  #power law
  x = ((x1^(n+1) - x0^(n+1))*y + x0^(n+1))^(1/(n+1))
  x_max = max(x)
  x = (-1*x+max(x))^log10(length(y))
  sign <-sample(c(-1,1),replace=T,length(x))
  effect_sizes <- x*sign
  pheno <- colSums(effect_sizes*t(genotype))
  range_main <- max(pheno)-min(pheno)
  range_factor <- range_main/range
  effect_sizes <- effect_sizes/range_factor
  filter<- sd(effect_sizes)
  effect_sizes[effect_sizes < sd(effect_sizes)] <- 0
  pheno <- colSums(effect_sizes*t(genotype))
  output<-(list(pheno,effect_sizes,filter))
  names(output)<-c("Phenotype","Effect_sizes","Filter")
  return(output)
}

enviroment_layer <- function(epi_layer,main_layer,n){
  genetic_effects<-epi_layer+main_layer
  baseline <- abs(min(genetic_effects)-max(genetic_effects))
  enviromental_layer <-rnorm(n=n,mean=0.1*baseline,sd=0.05*baseline)
  output <-list(baseline+enviromental_layer+genetic_effects,baseline,enviromental_layer)
  names(output) <- c("Phenotype","baseline","enviromental")
  return(output)
}








cor.prob <- function(X, dfr = nrow(X) - 2) {
  R <- cor(X)
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr / (1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
    cor.mat <- t(R)
  cor.mat[upper.tri(cor.mat)] <- NA
  cor.mat
}



library(wish)



epistasis_generation <- function(main_pheno,main_effect,genotype,n_snps,hera_exp){ 
  epi_pheno <- rep(0,length(main_pheno))
  total_effects<-sum(main_effect)/hera_exp
  likelyhood<-abs(main_effect %*% t(main_effect))
  likelyhood[likelyhood < min(main_effect)] <- 0
  likelyhood[lower.tri(likelyhood)] <-0 
  effect_vector <- sort(c(likelyhood[likelyhood > 0]),decreasing = T)
  index <- 1
  epi_total <- 0
  while(epi_total < total_effects & index < length(effect_vector)){
    epi_total <- epi_total+effect_vector[index]
    print(c(index,effect_vector[index],total_effects,epi_total))
    index <- index+1
  }
  likelyhood[likelyhood < effect_vector[index]] <- 0
  epi_effects<- likelyhood
  for (i in 1:n_snps){
    if (i+1 <= n_snps){
      for (j in (i+1):n_snps){
        epi_geno<-genotype[,i]*genotype[,j]      
        epi_pheno <- epi_pheno+epi_effects[i,j]*epi_geno
        
      }  
    }
  }
  output<-(list(epi_pheno,epi_effects,filter))
  names(output)<-c("Phenotype","Effect_sizes","Filter")
  return(output)
}


epi_pheno<-epistasis_layer(sampled_data,n_pheno=200)

phenotype<-enviroment_layer(epi_pheno[[1]],main_pheno[[1]],200)
hist(phenotype[[1]])

test_phenotype<-phenotype[[1]]
main_pheno[[1]]

sampled_data<-uniform_genotype_sampler(200,1000)

hist(main_pheno$Effect_sizes[main_pheno$Effect_sizes > 0])
hist(main_pheno$Phenotype)

sum(main_pheno$Effect_sizes > 0)
sum(main_pheno$Effect_sizes > sd(main_pheno$Effect_sizes))


sort(main_pheno$Effect_sizes)
sampled_data<-uniform_genotype_sampler(200,200)
epi_pheno <- epistasis_generation(main_pheno$Phenotype,sampled_data,200,0.1)
main_pheno<-main_layer(sampled_data,range = 20)
hist(main_pheno$Effect_sizes)
final_pheno<-enviroment_layer(epi_pheno$Phenotype,main_pheno$Phenotype,203)

length


library(WISH)

sampled_genotype<-genotype[,sample(651720,1000)]
sampled_genotype[is.na(sampled_genotype)]<-0
dim(sampled_genotype)

main_pheno<-main_layer(as.matrix(sampled_genotype),range = 20)
epi_pheno <- epistasis_generation(main_pheno$Phenotype,main_pheno$Effect_sizes,sampled_genotype,200,0.1)
hist(main_pheno$Effect_sizes)
hist(epi_pheno$Effect_sizes)
final_pheno<-enviroment_layer(epi_pheno$Phenotype,main_pheno$Phenotype,203)
test_corr<-epistatic.correlation(epi_pheno$Phenotype+main_pheno$Phenotype,sampled_genotype,threads = 20,test=F) 



test_corr<-epistatic.correlation(rnorm(203,100,50),sampled_genotype,threads = 20,test=F) 

adjusted<-(p.adjust(test_corr$Pvalue,method="bonferroni"))
sum(adjusted< 0.05, na.rm=T)


test_corr$Pvalues[is.na(test_corr$Pvalues)] <- 1
test_corr$Coefficients[is.na(test_corr$Coefficients)]<-0

test_corr$Coefficients[epi_pheno$Effect_sizes >0]

hist(test_corr$Pvalues[epi_pheno$Effect_sizes >0]) 

hist(epi_pheno$Effect_sizes[test_corr$Pvalues< 0.01]) 

cor(test_corr$Pvalues[epi_pheno$Effect_sizes >0],epi_pheno$Effect_sizes[epi_pheno$Effect_sizes >0])
     
cor(test_corr$Coefficients[epi_pheno$Effect_sizes >0],epi_pheno$Effect_sizes[epi_pheno$Effect_sizes >0])    


plot(test_corr$Coefficients[epi_pheno$Effect_sizes >0],epi_pheno$Effect_sizes[epi_pheno$Effect_sizes >0])    



dim(test_corr$Pvalues)  
dim(epi_pheno$Effect_sizes)  
sum(epi_pheno$Effect_sizes >0)  

cor(test_corr$Coefficients[epi_pheno$Effect_sizes >0],epi_pheno$Effect_sizes[epi_pheno$Effect_sizes >0])   

library(data.table) 
 
tped_variation<-fread("/data/nbg153/ADHD/WISH_files_newest/INDICES_1_2.tped",data.table = F) 
ped1<-fread("/data/nbg153/ADHD/WISH_files_newest/part1.ped",data.table = F) 
ped2<-fread("/data/nbg153/ADHD/WISH_files_newest/part2.ped",data.table = F) 
ped3<-fread("/data/nbg153/ADHD/WISH_files_newest/part3.ped",data.table = F) 
ped_variation<-cbind(ped1,ped2,ped3) 


genotype<-generate.genotype 

dim(genotype) 
  
sampled_genotype<-genotype[,sample(651720,200)] 

total_minors<-sum(tped_variation[,5:410] == 1) 
total_majors<-sum(tped_variation[,5:410] == 2) 
total_minors/(total_majors+total_minors) 

freq_minors<-rowSums(tped_variation[,5:410] == 1) 
freq_majors<-rowSums(tped_variation[,5:410] == 2) 

hist(freq_minors)
hist(freq_majors)
median(freq_minors)

