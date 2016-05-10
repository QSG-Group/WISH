#' Data format Import
#' @export
#' @param ped ?
#' @param tped true values
#' @return R-squared the fraction of explained sqaures of y by x
#' @details There is so much to be said
#' @examples
#' @import 
#'
data.import <- function(ped, tped,gwas,pvalue,id.select){
  a.s <- ped[,c(2,7:ncol(ped))]
  snp <- rep((as.character(tped[,2])),each=2)
  snp <- append(snp, "ID", after=0)
  colnames(a.s) <- c(snp)
  gwas.select <- as.data.frame(cbind(row.names(gwas),gwas@results$P1df))
  gwas.select$V2 <- as.numeric(as.character(gwas.select$V2))
  p.select <- gwas.select[gwas.select$V2 < pvalue,]
  snp.select <- p.select$V1
  snps<- rep(snp.select, each=2)
  snps <- append(as.character(snps), "ID", after=0)
  a.snps <- a.s[,snps]
  id.snps <- a.snps[a.snps$ID %in% id.select,]
  rownames(id.snps) <- id.snps[,1]
  id.snps <- id.snps[,2:ncol(id.snps)]
  return(id.snps)
}
genotype.data <- function(input){
  ind <- data.frame(matrix(c((1:ncol(input)),rep(NA, 2-ncol(input)%%2)),byrow=F,nrow=2))
  nonna <- ind[,sapply(ind, function(x) all(!is.na(x)))]
  genotype <- do.call(cbind,lapply(nonna,function(i)rowMeans(input[,i])))
  x <- colnames(input[1:(ncol(input)-1)])
  colnames(genotype) <- x[c(T,F)]
  genotype <- replace(genotype, genotype==0,"NA")
  genotype <- replace(genotype, genotype==2,3 )
  genotype <- replace(genotype, genotype==1.5,2 )
  genotype <- replace(genotype, genotype==1,1 )
  return(genotype)
}
create.matrix <- function(gwas,pvalue){
  gwas.select <- as.data.frame(cbind(row.names(gwas),gwas@results$P1df))
  gwas.select$V2 <- as.numeric(as.character(gwas.select$V2))
  p.select <- gwas.select[gwas.select$V2 < pvalue,]
  nsnps <- nrow(p.select)
  matrix <- matrix(NA,nrow=nsnps,ncol=nsnps)
  rownames(matrix) <- p.select$V1
  colnames(matrix) <- p.select$V1
  return(matrix)
}

epistatic.correlation <- function(gwas, pvalue, matrix, phenotype, genotype, trait){
  cand <- which(gwas@results$P1df <= pvalue)
  phenotype <- phenotype[with(phenotype,order(id)),]
  genotype <- genotype[rownames(genotype) %in% phenotype$id,]
  genotype <- genotype[,colnames(genotype) %in% rownames(gwas)]
    for (i in 1:(nrow(matrix)-1)){
    matrix[i, (i+1):nrow(matrix)] <- foreach(j=(i+1):nrow(matrix), .combine='rbind', .inorder=T, .verbose=F) %dopar% {
      phenotype$SNP1 = (genotype[,cand[i]]);
      phenotype$SNP2 = (genotype[,cand[j]]);
      correlation = asreml(trait~1+SNP1+SNP2+SNP1*SNP2,na.method.X="omit", na.method.y="omit",data=phenotype)
      return(correlation$coefficients$fixed[1])
    }
  }
}
