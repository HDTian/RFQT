\name{GetTree}
\alias{GetTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Fit and Get a Q Tree
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{GetTree} fits a data and construct a Q-tree. The fitted data is directly used to build tree (so no bootstrap). The tree information will be stored but not cleaned for further analysis. You may need \code{BootstrapTreeFitting} rather than this function.
}
\usage{
GetTree(dat,
        S = 5,
        Qthreshold = 3.84,
        rate = 1,
        SpecificM = NA,
        method = "DR",
        SoP = 10,
        howGX = "SpecificGX",
        Halve = FALSE,
        const = NA,
        endsize = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{dat}{a data frame. This data (without any bootstrap) will be fitted to construct a Q-tree.}
\item{S}{a postive integer indicates the maximum tree depth allowed. Default value is 5 (i.e. the data set is allowed to be splitted at most 5 times). S is also conneted to the argument endsize, both of which work in a similar way. One may only restrict one of them and leave another one be a flexible value (e.g. \code{S=100} or \code{endsize=0}). }
\item{Qthreshold}{A positive value used as the threshold value for the Q statistic. Default value is 3.0.}
\item{rate}{a value ranging from 0 to 1 indicates the proportion of the candidate covariates to be randomly cosidered in each split when construcing the Q-tree. When \code{rate=1}, it means all the candiate covariates will be considered.}
\item{SpecificM}{a vector indicates the index of candidate variables that can be considered in spliiting. Tbe default value is \code{NA}, which means all candidate variables are possible to be considered.}
\item{method}{a character indicates the stratification method used for constructing a Q-tree. There are currently three methods: the Doubly-ranked method (\code{method='DR'}), the residual methods (\code{method='Residual'}) and the naive method (recognized by any character except \code{'DR'} and \code{'Residual'}). The default method is \code{'DR'}. }
\item{SoP}{a postive integer (>=2) indicates the size of pre-straum. Only applicable for the doubly-ranked stratification. Default value is 10. }
\item{howGX}{a character indicates the way to calculate the instrument-exposure associations. Two ways are currenly allowed: the instrument-exposure associations estimated seperately at each stratum (\code{howGX='SpecificGX'}), or use a fixed instrument-exposure association (\code{howGX='const')}}
\item{Halve}{Logical. Indictes whether to split the node by half:half (i.e. 5:5). Default value is \code{FALSE}.}
\item{const}{a value indicates the fixed instrument-exposure association value used for constructing the Q-tree. Only applicable for \code{howGX='SpecificGX'}. }
\item{endsize}{a positive integer value indicates the minimual size of the node of Q-tree allowed. One can also control the endsize via \code{S}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{dat} must be a data frame that can be regonized by the functions. The first four columns must be the individual IDs, the (one-dimensional) instrument, the exposure (should be \code{NA} if not available), and the outcome, respectively. The (high-dimensional) covariates are then follwing. Note that all the functions may be working with incorrect column ordering, but the results are misleading. DO make sure the column variable order is correct.

If simulated data is used, there should be a end column named as true_STE. That is, \code{dat$true_STE} is not \code{NULL}.
}
\value{
\code{GetTree} returns a data set with the same row dimension as the inputed data \code{dat}. The returned data is 'augmented' with tree information including the tree end node index, and the splitting information in each split.
}
\references{
\cite{Burgess, S., Davies, N. M., & Thompson, S. G. (2014). "Instrumental variable analysis with a nonlinear exposure–outcome relationship". Epidemiology (Cambridge, Mass.), 25(6), 877. (Residual stratification)}

\cite{Tian, H., Mason, A. M., Liu, C., & Burgess, S. (2023). "Relaxing parametric assumptions for non-linear Mendelian randomization using a doubly-ranked stratification method". PLOS Genetics, 19(6), e1010823. (Doubly-ranked stratification)}
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

res<-getDat()
odat<-res$traning.set  #training set

tree_result<-GetTree(odat)
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
