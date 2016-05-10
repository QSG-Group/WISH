\name{create.matrix}
\alias{create.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
create.matrix(gwas, pvalue)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gwas}{
%%     ~~Describe \code{gwas} here~~
}
  \item{pvalue}{
%%     ~~Describe \code{pvalue} here~~
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
function (gwas, pvalue) 
{
    gwas.select <- as.data.frame(cbind(row.names(gwas), gwas@results$P1df))
    gwas.select$V2 <- as.numeric(as.character(gwas.select$V2))
    p.select <- gwas.select[gwas.select$V2 < pvalue, ]
    nsnps <- nrow(p.select)
    matrix <- matrix(NA, nrow = nsnps, ncol = nsnps)
    rownames(matrix) <- p.select$V1
    colnames(matrix) <- p.select$V1
    return(matrix)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
