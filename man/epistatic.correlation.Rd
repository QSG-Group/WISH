\name{epistatic.correlation}
\alias{epistatic.correlation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
epistatic.correlation(gwas_id, gwas_p, pvalue = 1e-04, phenotype, genotype, parallel = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gwas_id}{
%%     ~~Describe \code{gwas_id} here~~
}
  \item{gwas_p}{
%%     ~~Describe \code{gwas_p} here~~
}
  \item{pvalue}{
%%     ~~Describe \code{pvalue} here~~
}
  \item{phenotype}{
%%     ~~Describe \code{phenotype} here~~
}
  \item{genotype}{
%%     ~~Describe \code{genotype} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
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
function (gwas_id, gwas_p, pvalue = 1e-04, phenotype, genotype, 
    parallel = 20) 
{
    registerDoParallel(parallel)
    gwas <- cbind(gwas_id, gwas_p)
    gwas <- as.data.frame(gwas)
    gwas$gwas_p <- as.numeric(as.character(gwas$gwas_p))
    gwas <- gwas[gwas$gwas_id \%in\% colnames(genotype), ]
    cand <- gwas[gwas$gwas_p < pvalue, ]
    snp_matrix <- matrix(NA, nrow = nrow(cand), ncol = nrow(cand))
    rownames(snp_matrix) <- cand$gwas_id
    colnames(snp_matrix) <- cand$gwas_id
    cand <- which(gwas$gwas_p <= pvalue)
    for (i in 1:(nrow(snp_matrix) - 1)) {
        snp_matrix[i, (i + 1):nrow(snp_matrix)] <- foreach(j = (i + 
            1):nrow(snp_matrix), .combine = "rbind", .inorder = T, 
            .verbose = F) \%dopar\% {
            phenotype$SNP1 = as.numeric((genotype[, cand[i]]))
            phenotype$SNP2 = as.numeric((genotype[, cand[j]]))
            tmp_correlation = asreml(OI ~ 1 + SNP1 + SNP2 + SNP1 * 
                SNP2, na.method.X = "omit", na.method.Y = "omit", 
                data = phenotype)
            return(tmp_correlation$coefficients$fixed[1])
        }
    }
    new <- t(snp_matrix)
    diag(new) <- 1
    new[upper.tri(new)] <- snp_matrix[upper.tri(snp_matrix)]
    return(new)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
