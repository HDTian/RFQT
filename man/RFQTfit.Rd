\name{RFQTfit}
\alias{RFQTfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a Random Forest of Q Trees (RFQT)
}
\description{
\code{RFQTfit} fits a data set (the training data) by the random forest of Q trees (RFQT), and return the summary results.
}
\usage{
RFQTfit(odat,
        vdat = NA,
        Nb = 1,
        S = 5,
        honest = FALSE,
        rate = 0.4,
        SingleM = FALSE,
        Qthreshold = 3.84,
        method = "DR",
        SoP = 10,
        howGX = "SpecificGX",
        Halve = FALSE,
        endsize = 1000,
        const = NA,
        Cores = NA,
        trackfile=NA
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
  \item{Nb}{
a positive integer indicates the number of bootstrap (i.e. the number of Q-trees)
}
  \item{S}{
a postive integer indicates the maximum tree depth allowed. Default value is 5 (i.e. the data set is allowed to be splitted at most 5 times). \code{S} is also conneted to the argument \code{endsize}, both of which work in a similar way. One may only restrict one of them and leave another one be a flexible value (e.g. \code{S=100} or \code{endsize=0}).
}
\item{honest}{
logical value indicates whether to use the honest estimation style or not. Default value is \code{FALSE}. When honest estimation is used, the tree construction and the leaf-specific estimates will be obtained in two seperate samples of the (boostrapped) training data.
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
\item{Cores}{
a positive integer indicates the number of cores should be used when running random forest in parallel. Default value is the number of CPU cores minus one. You can check the number of CPU cores via \code{detectCores()}.
}
\item{trackfile}{
the path string indicates the file to store the running information in parallel. 
}

}
\details{
\code{RFQTfit} is an integrated function combining other functions, including \code{GetTree}, \code{BootstrapTreeFitting}, \code{getPredict}, \code{getVI}, etc. The fitting automatically uses parallel computation.

The training data \code{Odat} and the testing data \code{Vdat} should be of the same form and be a data frame that can be regonized by the functions. The first four columns must be the individual IDs, the (one-dimensional) instrument, the exposure (should be \code{NA} if not available), and the outcome, respectively. The (high-dimensional) covariates are then follwing. Note that all the functions may be working with incorrect column ordering, but the results are misleading. DO make sure the column variable order is correct.

}
\value{
\code{RFQTfit} shows all the parameter information used for random forest fitting in the begining, and returns in the end a list consisting of the following components:
\item{RES}{
The direct fitting results of each Q-tree, which contains \code{end_node_information}, \code{OOB_predict}, \code{v_predict}, \code{vi1}, \code{vi2}, \code{ts1} and \code{ts2} (same as the results of \code{BootstrapTreeFitting}) for each Q-tree. Note that \code{RES} is multi-dimensional. If one wishes to obtain the individual predicted effects of the testing set according to i-th Q-tree, should run \code{RES[3,i]$v_predict}.
}
\item{MSE_OOB}{
the \code{Nb}-length vector of MSE results for OOB samples. If one determine the number of Q-trees to be \code{Nb}, simply use the final vlue of \code{MSE_OOB} as the MSE.
}
\item{MSE_test}{
the \code{Nb}-length vector of MSE results for test samples (if the testing set exists). If one determine the number of Q-trees to be \code{Nb}, simply use the final vlue of \code{MSE_test} as the MSE.
}
\item{Predict_OOB}{
the matrix of predicted effects for OOB samples, where the row corresponds to the sample size (the testing data) and the column coresponds to the number of Q-trees. Each column represents the individual predicted effects with the present number of Q trees. If one determine the number of Q-trees to be \code{Nb}, simply use the final colunm of \code{Predict_OOB} as the individual predicted effect vector. 
}
\item{Predict_test}{
the matrix of predicted effects for testing set, where the row corresponds to the sample size (the testing data) and the column coresponds to the number of Q-trees. Each column represents the individual predicted effects with the present number of Q trees. If one determine the number of Q-trees to be \code{Nb}, simply use the final colunm of \code{Predict_test} as the individual predicted effect vector. 
}
\item{VI1}{
the vector indicates the type I (i.e. with labels) variable importacne (VI) measurements for all the candidate covariates (with the same ordering of the candidate covariates to the original fitting data).
}
\item{VI2}{
the vector indicates the type II (i.e. without labels) variable importacne (VI) measurements for all the candidate covariates (with the same ordering of the candidate covariates to the original fitting data).
}
\item{TS}{
the vector of the permutation test statistic values. The two values represent the values of TS1 and TS2, respectively.
}

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
#res<-getDat() #simulated data  #the deflaut setting: scenario='A' and SoM=0.5
#odat<-res$traning.set  
#vdat<-res$testing.set  

#ALLRES<-RFQTfit(odat)
 
#ALLRES<-RFQTfit(odat,vdat,Nb=200)


}

