\name{GetIndex}
\alias{GetIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get the Best Covariate Index
}
\description{
\code{GetIndex} gives the best covariate index according to the inputed data set. The best covariate is one out of the candidate covariates considered that gives the largest Q statistic value.
}
\usage{
GetIndex(dat_current,
          JJ,
          rate = 1,
          SpecificM = NA,
          method = "DR",
          SoP = 10,
          howGX = "SpecificGX",
          Halve = FALSE,
          const = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat_current}{ A data frame. This data frame has the first four columns: individual IDs, the instrument, the expoure, and the outcome with the candidate covariates following.
}
  \item{JJ}{A positive integer indicates the total number of the candiate covariates.
}
  \item{rate}{a value ranging from 0 to 1 indicates the proportion of the candidate covariates to be randomly cosidered in each split when construcing the Q-tree.}
  \item{SpecificM}{a vector indicates the index of candidate variables that can be considered in spliiting. Tbe default value is \code{NA}, which means all candidate variables are possible to be considered.}
  \item{method}{a character indicates the stratification method used for constructing a Q-tree. There are currently three methods: the Doubly-ranked method (\code{method='DR'}), the residual methods (\code{method='Residual'}) and the naive method (recognized by any character except \code{'DR'} and \code{'Residual'}). The default method is \code{'DR'}.
}
  \item{SoP}{a postive integer (>=2) indicates the size of pre-straum. Only applicable for the doubly-ranked stratification. Default value is 10.
}
  \item{howGX}{a character indicates the way to calculate the instrument-exposure associations. Two ways are currenly allowed: the instrument-exposure associations estimated seperately at each stratum (\code{howGX='SpecificGX'}), or use a fixed instrument-exposure association (\code{howGX='const')}
}
\item{Halve}{Logical. Indictes whether to split the node by half:half (i.e. 5:5). Default value is \code{FALSE}.}
  \item{const}{a value indicates the fixed instrument-exposure association value used for constructing the Q-tree. Only applicable for \code{howGX='SpecificGX'}.
}
}
\details{
\code{GetIndex} helps to decide the splitting covariate at the present node when constructing a Q-tree.
}
\value{
\code{GetIndex} returns a vector with the following three values:'Candidate.index', 'Q.value' , 'node.size'
\item{Candidate.index }{The index (position rank) of the chosen covariate.}
\item{Q.value}{The corresponding Q statistic value of the chosen covariate.}
\item{minimal.node.size.after.split}{The mimimal sample size of the sun-nodes after splitting.}
\item{split.style}{Three values denoted by 1,2,3 representing the splitting style used for the final covariate. They are 3:7, 5:5, 7:3, respectively. }  
}

\references{
\cite{Burgess, S., Davies, N. M., & Thompson, S. G. (2014). "Instrumental variable analysis with a nonlinear exposure–outcome relationship". Epidemiology (Cambridge, Mass.), 25(6), 877. (Residual stratification)}

\cite{Tian, H., Mason, A. M., Liu, C., & Burgess, S. (2022). "Relaxing parametric assumptions for non-linear Mendelian randomization using a doubly-ranked stratification method". bioRxiv, 2022-06. (Doubly-ranked stratification)}
}
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
library(MendelianRandomization)
library(tidyverse) 
library(data.table) 

set.seed(60)
res<-getDat() #simulated data #the deflaut setting: scenario='A' and SoM=0.5
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set
GetIndex(odat,JJ=20)
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
