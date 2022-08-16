\name{doKernel}
\alias{doKernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
doKernel
}
\description{
Obtain the kernel gene expression profile given gene expression matrix, pseudotime and cell weights.
}
\usage{
dokernel(count,pseudotime, cellWeights, nGenes, nKnots, nPoints, fdr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{count}{
    an expression data frame, each column is a             query cell.}
    \item{pseudotime}{
    A vector or dataframe of the cell pseudotime. We suggest to use the output of psuedotime analysis methods such as Slingshot, monocle, PAGA, etc.}
    \item{cellWeights}{
    A vector or dataframe implying the branch of the cell.}
    \item{nGenes}{
    In order to save computation time, one can define number of top variable genes to be considered for the candidate genes of kernel expresion profile.}
    \item{nKnots}{
    See tradeSeq for optimizing the parameter nKnots for inferring gene expression dynamics associated with pseudotime.}
     \item{nPoints}{
    Number of psudo timepoints in the kernel profile. Default value is 100}
    \item{fdr}{
    Adjusted P value cutoff for selecting kernel genes. Default value is 0.05}
 }
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A kernel expression profile depicting the most important dynamic gene expressions along the trajectory.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
\author{
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
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
