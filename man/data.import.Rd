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
data.import(ped, tped, gwas, pvalue, id.select)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{
%%     ~~Describe \code{ped} here~~
}
  \item{tped}{
%%     ~~Describe \code{tped} here~~
}
  \item{gwas}{
%%     ~~Describe \code{gwas} here~~
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
function (ped, tped, gwas, pvalue, id.select) 
{
    a.s <- ped[, c(2, 7:ncol(ped))]
    snp <- rep((as.character(tped[, 2])), each = 2)
    snp <- append(snp, "ID", after = 0)
    colnames(a.s) <- c(snp)
    gwas.select <- as.data.frame(cbind(row.names(gwas), gwas@results$P1df))
    gwas.select$V2 <- as.numeric(as.character(gwas.select$V2))
    p.select <- gwas.select[gwas.select$V2 < pvalue, ]
    snp.select <- p.select$V1
    snps <- rep(snp.select, each = 2)
    snps <- append(as.character(snps), "ID", after = 0)
    a.snps <- a.s[, snps]
    id.snps <- a.snps[a.snps$ID \%in\% id.select, ]
    rownames(id.snps) <- id.snps[, 1]
    id.snps <- id.snps[, 2:ncol(id.snps)]
    return(id.snps)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
