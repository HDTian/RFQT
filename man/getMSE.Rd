\name{getMSE}
\alias{getMSE}
\title{
Get MSE Values for the OOB Data or the Testing Data with Fitted Q-Tree or Random Forest of Q-Trees (RFQT)
}
\description{
\code{getMSE} calculates the MSE value for RFQT (or a Q-tree if RFQT only contains one single tree)
}
\usage{
getMSE(RES, 
       indicator = 1
       )
}
\arguments{
  \item{RES}{
a fitting list result from \code{BootstrapTreeFitting} or \code{RFQTfit(odat,vdat)$RES}.
}
  \item{indicator}{
numeric indicator to judge which type of MSE should be calculated. \code{indicator=1} for Out-of-Bag (OOB) data and \code{indicator=2} for th testing data.  Defalut value is \code{indicator=1}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\code{getMSE} returns a vector of MSE with the increase of the number of Q-trees. The stable MSE values should indicate an appropriate number of Q-trees.
}
\references{
%% ~put references to the literature/web site here ~
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
set.seed(60)
res<-getDat() #simulated data  #the deflaut setting: scenario='A' and SoM=0.5
res<-getDat(label=FALSE) 
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set
 
#When running RFQT with mutiple Q trees/bootstrap - use parallel computation
Nb<-7 #how many trees in the forest? e.g. 5 trees
cl<-makeCluster(detectCores()-1)
clusterEvalQ(cl=cl , expr=library(dplyr)) 
clusterEvalQ(cl=cl , expr=library(MendelianRandomization) )
clusterExport(  cl=cl ,  varlist=c('odat','vdat',
                                    'GetTree', 'GetNindex', 'GetIndex' )  )
RES<-parSapply(   cl ,  1:Nb, BootstrapTreeFitting  )
stopCluster(cl)
dim(RES)#7 Nb
 
 
getMSE(RES,indicator=1 )
getMSE(RES,indicator=2 )
}

