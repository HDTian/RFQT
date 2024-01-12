\name{getVI}
\alias{getVI}
\title{
Get Varibale Importance Measurements with Fitted Q-Tree or Random Forest of Q-Trees (RFQT)
}
\description{
\code{getVI} calculates the Varibale Importance (VI) measurements for a random forest of Q-trees (RFQT).
}
\usage{
getVI(RES, 
      VItype = 2
      )
}
\arguments{
  \item{RES}{
a fitting list result from \code{BootstrapTreeFitting} or \code{RFQTfit(odat,vdat)$RES}.
}
  \item{VItype}{
numeric indicator to judge which type of VI measurement should be calculated. \code{VItype=1} for the classical way with the known individual treatment effects (i.e. label) and \code{VItype=2} for case that the individual treatment effects are unavailable (as in most real application).  Defalut value is \code{VItype=2}.
}
}
\details{
See the VI measurement part of the original paper for more detailed VI introduction and the algorithm.
}
\value{
\code{getVI} retuns a vector indicating the variable importacne (VI) measurements for all the candidate covariates (with the same ordering of the candidate covariates to the original fitting data).
}

\examples{
library(MendelianRandomization)
library(tidyverse) 
library(data.table) 
library(parallel)

set.seed(60)
res<-getDat() #simulated data  #the deflaut setting: scenario='A' and SoM=0.5
res<-getDat(label=FALSE) 
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set
 
##When running RFQT with mutiple Q trees/bootstrap - use parallel computation
Nb<-7 #how many trees in the forest? e.g. 5 trees
cl<-makeCluster(2)
clusterEvalQ(cl=cl , expr=library(dplyr)) 
clusterEvalQ(cl=cl , expr=library(MendelianRandomization) )
clusterExport(  cl=cl ,  varlist=c('odat','vdat',
                                    'GetTree', 'GetNindex', 'GetIndex' )  )
RES<-parSapply(   cl ,  1:Nb, BootstrapTreeFitting  )
stopCluster(cl)
dim(RES)#7 Nb

getVI(RES,VItype=1 )
getVI(RES,VItype=2 )
}

