\name{DTfit}
\alias{DTfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a Classical Single Decision Tree (DT)
}
\description{
\code{DTfit} fits a data set (the training data) by the classical Q tree, and return the summary results. \code{DTfit} works similarly as \code{GetTree} and \code{BootstrapTreeFitting}, but not designed or used for future random forest fitting. If you wish to fit data via single tree rather than random forest, use \code{DTfit}. 
}
\usage{
RFQTfit(Odat=odat,  
        Vdat=NA,
        honest=FALSE,
        S=5, 
        rate=1, 
        SingleM=FALSE, 
        Qthreshold=3.84, 
        method='DR',
        SoP=10, 
        howGX='SpecificGX',
        Halve=FALSE, 
        endsize=1000,
        const=NA)
        )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{odat}{
a data frame. This is recognized as the trainging data set.
}
  \item{vdat}{
a data frame. This is recognized as the testing or validation data set. The default value is \code{NA}, indicating that no testing data will be considered. 
}
\item{honest}{
logical value indicates whether to use the honest estimation style or not. Default value is \code{FALSE}. When honest estimation is used, the tree construction and the leaf-specific estimates will be obtained in two seperate samples of the (boostrapped) training data.
}
  \item{S}{
a postive integer indicates the maximum tree depth allowed. Default value is 5 (i.e. the data set is allowed to be splitted at most 5 times). \code{S} is also conneted to the argument \code{endsize}, both of which work in a similar way. One may only restrict one of them and leave another one be a flexible value (e.g. \code{S=100} or \code{endsize=0}).
}
  \item{rate}{
a value ranging from 0 to 1 indicates the proportion of the candidate covariates to be randomly cosidered in each split when construcing the Q-tree.
}
\item{SingleM}{
logical value indicates whether to use one single covariate variable at all times. Default value is \code{FALSE}. If \code{SingleM=TRUE}, the covariate used are the best variable indicated by \code{GetIndex}.
}
  \item{Qthreshold}{
a positive value used as the threshold value for the Q statistic. Default value is 3.84. One can set \code{Qthreshold=0} to ignore the Q-related stopping rule.
}
  \item{method}{
a character indicates the stratification method used for constructing a Q-tree. There are currently three methods: the Doubly-ranked method (\code{method='DR'}), the residual methods (\code{method='Residual'}) and the naive method (recognized by any character except \code{'DR'} and \code{'Residual'}). The default method is \code{'DR'}.
}
  \item{SoP}{
a postive integer (>=2) indicates the size of pre-straum. Only applicable for the doubly-ranked stratification. Default value is 10.
}
  \item{howGX}{
a character indicates the way to calculate the instrument-exposure associations. Two ways are currenly allowed: the instrument-exposure associations estimated seperately at each stratum (\code{howGX='SpecificGX'}), or use a fixed instrument-exposure association (\code{howGX='const'}).
}
\item{Halve}{
logical. Indictes whether to split the node by half:half (i.e. 5:5). Default value is \code{FALSE}.
}
  \item{endsize}{
a positive integer value indicates the minimual size of the node of Q-tree. \code{S} and \code{endsize} work in a similar way.
}
  \item{const}{
a value indicates the fixed instrument-exposure association value used for constructing the Q-tree. Only applicable for \code{howGX='SpecificGX'}.
}


}
\details{
Like \code{RFQTfit}, \code{DTfit} is an integrated function combining other functions, including \code{getPredict} and \code{getVI}.

The training data \code{Odat} and the testing data \code{Vdat} should be of the same form and be a data frame that can be regonized by the functions. The first four columns must be the individual IDs, the (one-dimensional) instrument, the exposure (should be \code{NA} if not available), and the outcome, respectively. The (high-dimensional) covariates are then follwing. Note that all the functions may be working with incorrect column ordering, but the results are misleading. DO make sure the column variable order is correct.

}
\value{
\code{DTfit} returns a list consisting of the following components:
\item{end_node_information}{
the end nodes (leaves) information of the fitted Q tree, which contains the end-node index, the MR estimates, and the sample size proportion.
}
\item{v_predicted}{the predicted effects of the testing samples accrding to the Q-tree.}
\item{ts1}{the value of the permutation test with statistic S1.}
\item{ts2}{the value of the permutation test with statistic S2.}
\item{MSE}{the MSE results for the testing samples (if both of the testing set and the label exist).}
\item{rdat}{the tree fitted data.}
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
#set.seed(60)
#res<-getDat() #simulated data  #the default setting: scenario='A' and SoM=0.5
#odat<-res$traning.set  #training set
#vdat<-res$testing.set  #testing set
 
#DTRES<-DTfit(odat,vdat)
 
#DTRES
 
##get (testing set) MSE manually
#mean( ( DTRES$v_predict - vdat$true_STE   )^2   )


}

