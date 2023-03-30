\name{GetNindex}
\alias{GetNindex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get the Q Tree End Node Index for choosen Samples
}
\description{
\code{GetNindex} helps to match the Q-tree end node for any inputed samples, and therefore enables to make further analysis like effect prediction for these inputed samples.
}
\usage{
GetNindex(M, rdat, S = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A matrix or data frame contains the candidate covariate information.
}
  \item{rdat}{
A data frame. It is the Q-tree informaton result and obtained by \code{rdat<-GetTree()}.
}
  \item{S}{A positive integer indicates the max depth of the tree (consistent with \code{rdat}).
%%     ~~Describe \code{S} here~~
}
}
\details{
The inputed data set \code{M} must contain the same candidate covariate information (with the same order) as the training data for the Q-tree. In addition, the first four columns (individual IDs, the instrument, the exposure, the outcome) need to be removed, so that the fist column starts from the first candidate covariate.
}
\value{
\code{GetNindex} returns a vector with the same length of the row number of the inputed data \code{M}. Each element represents the predicted effect according to the reference Q-tree (information expressed by \code{rdat}).
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Haodong Tian
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(MendelianRandomization)
library(tidyverse)
library(data.table)

res<-getDat()
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set

rdat<-GetTree(odat)

GetNindex(  odat[,-(1:4)] ,rdat )
}
