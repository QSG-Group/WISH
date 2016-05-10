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
epistatic.correlation(gwas, pvalue, matrix, phenotype, genotype, trait)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gwas}{
%%     ~~Describe \code{gwas} here~~
}
  \item{pvalue}{
%%     ~~Describe \code{pvalue} here~~
}
  \item{matrix}{
%%     ~~Describe \code{matrix} here~~
}
  \item{phenotype}{
%%     ~~Describe \code{phenotype} here~~
}
  \item{genotype}{
%%     ~~Describe \code{genotype} here~~
}
  \item{trait}{
%%     ~~Describe \code{trait} here~~
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
function (gwas, pvalue, matrix, phenotype, genotype, trait) 
{
    cand <- which(gwas@results$P1df <= pvalue)
    phenotype <- phenotype[with(phenotype, order(id)), ]
    genotype <- genotype[rownames(genotype) \%in\% phenotype$id, 
        ]
    genotype <- genotype[, colnames(genotype) \%in\% rownames(gwas)]
    for (i in 1:(nrow(matrix) - 1)) {
        matrix[i, (i + 1):nrow(matrix)] <- foreach(j = (i + 1):nrow(matrix), 
            .combine = "rbind", .inorder = T, .verbose = F) \%dopar\% 
            {
                phenotype$SNP1 = (genotype[, cand[i]])
                phenotype$SNP2 = (genotype[, cand[j]])
                correlation = asreml(trait ~ 1 + SNP1 + SNP2 + 
                  SNP1 * SNP2, na.method.X = "omit", na.method.y = "omit", 
                  data = phenotype)
                return(correlation$coefficients$fixed[1])
            }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
