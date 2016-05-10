\name{genotype.data}
\alias{genotype.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
genotype.data(input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
%%     ~~Describe \code{input} here~~
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
function (input) 
{
    ind <- data.frame(matrix(c((1:ncol(input)), rep(NA, 2 - ncol(input)\%\%2)), 
        byrow = F, nrow = 2))
    nonna <- ind[, sapply(ind, function(x) all(!is.na(x)))]
    genotype <- do.call(cbind, lapply(nonna, function(i) rowMeans(input[, 
        i])))
    x <- colnames(input[1:(ncol(input) - 1)])
    colnames(genotype) <- x[c(T, F)]
    genotype <- replace(genotype, genotype == 0, "NA")
    genotype <- replace(genotype, genotype == 2, 3)
    genotype <- replace(genotype, genotype == 1.5, 2)
    genotype <- replace(genotype, genotype == 1, 1)
    return(genotype)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
