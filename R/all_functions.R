#library("GenABEL")
#library("asreml")
#library("doParallel")
#library("sna")

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
