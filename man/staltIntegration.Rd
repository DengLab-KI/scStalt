\name{staltIntegration}
\alias{staltIntegration}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
staltIntegration
}
\description{
Use the reference trajectory to integrate two datasets.
}
\usage{
staltIntegration(count1, count2, nPoints = 100, ref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{count1}{
Expression matrix of dataset 1.}
  \item{count2}{
Expression matrix of dataset 2
}
  \item{nPoints}{
The pseudotime partritions to use for aligning different datasets
}
  \item{ref}{
Expression matrix of the reference trajectory
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The corrected counts matrix of the combined dataset 1 and 2
}

\author{
Qing Luo
}




