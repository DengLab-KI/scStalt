\name{stageAssign}
\alias{stageAssign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
stageAssign
}
\description{
Assign query cells to the stages in the kernel
}
\usage{
stageAssign(count,kernel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{count}{
    an expression data frame, each column is a             query cell and the genes are in rows in the           same sequence as in the kernel.}
    \item{kernel}{
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
