\name{BootstrapTreeFitting}
\alias{BootstrapTreeFitting}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
One Bootstrap Q Tree Fitting
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{BootstrapTreeFitting} fits the Q-tree for a bootstrap data of the training set (with user-specific seed) and store the cleaned results. If you do not need results but the original (dirty) Q-tree information, you may use \code{GetTree()}.
}
\usage{
BootstrapTreeFitting(seed=1,
                     Odat=odat,
                     Vdat=vdat,
                     honest=FALSE,
                     S=5,
                     rate=0.4,
                     SingleM=FALSE,
                     Qthreshold=3.84,
                     method='DR',
                     SoP=10,
                     howGX='SpecificGX',
                     Halve=FALSE,
                     endsize=1000,
                     const=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{seed}{seed value for reproducible and traceable results}
\item{Odat}{a data frame. This is recognized as the trainging data set.}
\item{Vdat}{a data frame. This is recognized as the testing or validation data set.}
\item{honest}{Logical value indicates whether to use the honest estimation style or not. Default value is \code{FALSE}. When honest estimation is used, the tree construction and the leaf-specific estimates will be obtained in two seperate samples of the (boostrapped) training data.}
\item{S}{a postive integer indicates the maximum tree depth allowed. Default value is 5 (i.e. the data set is allowed to be splitted at most 5 times). \code{S} is also conneted to the argument \code{endsize}, both of which work in a similar way. One may only restrict one of them and leave another one be a flexible value (e.g. \code{S=100} or \code{endsize=0}). }
\item{rate}{a value ranging from 0 to 1 indicates the proportion of the candidate covariates to be randomly cosidered in each split when construcing the Q-tree.}
\item{SingleM}{Logical value indicates whether to use one single covariate variable at all times. Default value is \code{FALSE}. If \code{SingleM=TRUE}, the covariate used are the best variable indicated by \code{GetIndex}.}
\item{Qthreshold}{A positive value used as the threshold value for the Q statistic. Default value is 3.0. One can set \code{Qthreshold=0} to ignore the Q-related stopping rule.}
\item{method}{a character indicates the stratification method used for constructing a Q-tree. There are currently three methods: the Doubly-ranked method (\code{method='DR'}), the residual methods (\code{method='Residual'}) and the naive method (recognized by any character except \code{'DR'} and \code{'Residual'}). The default method is \code{'DR'}. }
\item{SoP}{a postive integer (>=2) indicates the size of pre-straum. Only applicable for the doubly-ranked stratification. Default value is 10. }
\item{howGX}{a character indicates the way to calculate the instrument-exposure associations. Two ways are currenly allowed: the instrument-exposure associations estimated seperately at each stratum (\code{howGX='SpecificGX'}), or use a fixed instrument-exposure association (\code{howGX='const'})}
\item{Halve}{Logical. Indictes whether to split the node by half:half (i.e. 5:5). Default value is \code{FALSE}.}
\item{endsize}{a positive integer value indicates the minimual size of the node of Q-tree. \code{S} and \code{endsize} work in a similar way.}
\item{const}{a value indicates the fixed instrument-exposure association value used for constructing the Q-tree. Only applicable for \code{howGX='SpecificGX'}.}
%%     ~~Describe \code{x} here~~
}

\details{
%%  ~~ If necessary, more details than the description above ~~
The training data \code{Odat} and the testing data \code{Vdat} should be of the same form and be a data frame that can be regonized by the functions. The first four columns must be the individual IDs, the (one-dimensional) instrument, the exposure (should be \code{NA} if not available), and the outcome, respectively. The (high-dimensional) covariates are then follwing. Note that all the functions may be working with incorrect column ordering, but the results are misleading. Do make sure the column variable order is correct.

If simulated data is used, there should be a end column named as true_STE. That is, \code{odat$true_STE} is not \code{NULL}.

If the data does not contain the exposure information, the RFQT can still be fitted, where one should use \code{howGX='const'} with a external instrument-exposure association (usually from independent studies).
}
\value{
%%  ~Describe the value returned
\code{BootstrapTreeFitting} returns a list with the following components:

\item{end_node_information}{the end nodes (leaves) information of the fitted Q tree, which contains the end-node index, the MR estimates, and the sample size proportion.}
\item{OOB_predict}{the predicted effects of the OOB samples accrding to the Q-tree. Values will be 0 if the individual is not out of bag. }
\item{v_predict}{the predicted effects of the testing samples accrding to the Q-tree.}
\item{vi1}{the vector of the variable importance (VI) measurements for all the candidate covariates considered.}
\item{vi2}{the vector of the variable importance (VI) measurements for all the candidate covariates considered. This VI measurements do not need the individual true effect information, therefore suitable for real application data.}
\item{ts1}{the value of the permutation test with statistic S1.}
\item{ts2}{the value of the permutation test with statistic S2.}

When the real data is fitted, \code{vi1} is not applicable. When only traning data is available, \code{v_predict} is not applicable.
%% ...
}
\references{
\cite{Burgess, S., Davies, N. M., & Thompson, S. G. (2014). "Instrumental variable analysis with a nonlinear exposure–outcome relationship". Epidemiology (Cambridge, Mass.), 25(6), 877. (Residual stratification)}

\cite{Tian, H., Mason, A. M., Liu, C., & Burgess, S. (2023). "Relaxing parametric assumptions for non-linear Mendelian randomization using a doubly-ranked stratification method". PLOS Genetics, 19(6), e1010823. (Doubly-ranked stratification)}
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
library(MendelianRandomization)
library(tidyverse) 
library(data.table) 
library(parallel)

set.seed(60)
res<-getDat() #simulated data  #the deflaut setting: scenario='A' and SoM=0.5
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set


#When running RFQT with mutiple Q trees/bootstrap - use parallel computation
Nb<-5 #how many trees in the forest? e.g. 5 trees
cl<-makeCluster(2)#规定用多少并行clusters/线程
clusterEvalQ(cl=cl , expr=library(dplyr))  
clusterEvalQ(cl=cl , expr=library(MendelianRandomization) )
clusterExport(  cl=cl ,  varlist=c( 'odat', 'vdat',
                                    'GetTree', 'GetNindex', 'GetIndex' )  )#or any other arguments
RES<-parSapply(   cl ,  1:Nb, BootstrapTreeFitting  ) 
stopCluster(cl)
dim(RES)#7 Nb
##If you wish to use your own parameters rather than the default parameters, try:
#general exmaple
user_BootstrapTreeFitting<-function(seed){
  RES<-BootstrapTreeFitting(seed,
                            honest=my.honest,
                            S=my.S,
                            JJ=my.JJ,   #or any partial of the arguments
                            rate=my.rate,
                            Qthreshold=my.Qthreshold,
                            method=my.method,
                            SoP=my.SoP,
                            howGX=my.howGX,
                            endsize=my.endsize)
  return(RES)
}
#specific exmaple
user_BootstrapTreeFitting<-function(seed){
  RES<-BootstrapTreeFitting(seed,SoP=20)
  return(RES)
}
###parallel computation  -> RFQT
Nb<-5 #how many trees in the forest?
cl<-makeCluster(2)
clusterEvalQ(cl=cl , expr=library(dplyr)) 
clusterEvalQ(cl=cl , expr=library(MendelianRandomization) )
clusterExport(  cl=cl ,  varlist=c( 'odat', 'vdat',
                                    'GetTree', 'GetNindex', 'GetIndex' , 'BootstrapTreeFitting')  )
RES<-parSapply(   cl ,  1:Nb, user_BootstrapTreeFitting  )
stopCluster(cl)

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
