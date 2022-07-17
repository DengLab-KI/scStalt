\name{stageAssign}
\alias{stageAssign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
stageAssign
}
\description{
Assign query cells to the stages in the kernal
}
\usage{
stageAssign(count,kernal)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{count}{
    an expression data frame, each column is a             query cell and the genes are in rows in the           same sequence as in the kernal.}
    \item{kernal}{
    an expression data frame, each column is a            smoothed time point and each row is one gene.     }
 }
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame with the assigned stage and score for each query cell.
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
