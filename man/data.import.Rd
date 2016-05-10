\name{data.import}
\alias{data.import}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
data.import(ped, tped, gwas_id, gwas_p, pvalue = 0.05, id.select = ped[, 2])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{
%%     ~~Describe \code{ped} here~~
}
  \item{tped}{
%%     ~~Describe \code{tped} here~~
}
  \item{gwas_id}{
%%     ~~Describe \code{gwas_id} here~~
}
  \item{gwas_p}{
%%     ~~Describe \code{gwas_p} here~~
}
  \item{pvalue}{
%%     ~~Describe \code{pvalue} here~~
}
  \item{id.select}{
%%     ~~Describe \code{id.select} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ped, tped, gwas_id, gwas_p, pvalue = 0.05, id.select = ped[, 
    2]) 
{
    a.s <- ped[, c(2, 7:ncol(ped))]
    snp <- rep((as.character(tped[, 2])), each = 2)
    snp <- append(snp, "ID", after = 0)
    colnames(a.s) <- c(snp)
    gwas.select <- cbind(gwas_id, gwas_p)
    gwas.select <- as.data.frame(gwas.select)
    gwas.select$gwas_p <- as.numeric(as.character(gwas.select$gwas_p))
    p.select <- gwas.select[gwas.select$gwas_p < pvalue, ]
    snp.select <- p.select$gwas_id
    snps <- rep(snp.select, each = 2)
    snps <- append(as.character(snps), "ID", after = 0)
    a.snps <- a.s[, snps]
    rownames(a.snps) <- a.snps[, 1]
    input <- a.snps[, 2:ncol(a.snps)]
    ind <- data.frame(matrix(c((1:ncol(input)), rep(NA, 2 - ncol(input)\%\%2)), 
        byrow = F, nrow = 2))
    nonna <- ind[, sapply(ind, function(x) all(!is.na(x)))]
    genotype <- do.call(cbind, lapply(nonna, function(i) rowMeans(input[, 
        i])))
    x <- colnames(input[1:(ncol(input) - 1)])
    colnames(genotype) <- x[c(T, F)]
    genotype <- replace(genotype, genotype == 0, "NA")
    return(as.data.frame(genotype))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
