\name{getPredict}
\alias{getPredict}

\title{
Get predict effects for the OOB Data or the Testing Data with Fitted Random Forest of Q-Trees (RFQT)
}
\description{
\code{getPredict} calculates the predicted individual treatment effect estimates for RFQT, either for the OOB data or the testing data.
}
\usage{
getPredict(RES, 
           indicator = 1
           )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RES}{
a fitting list result from \code{BootstrapTreeFitting} or \code{RFQTfit(odat,vdat)$RES}.
}
  \item{indicator}{
numeric indicator to judge which type of predicted effect should be presented. \code{indicator=1} for Out-of-Bag (OOB) data and \code{indicator=2} for th testing data. Defalut value is \code{indicator=1}.
}
}
\details{
See the relevant part of the original paper for more details of how to predict an individual effect according to a fitted Q tree (therefore a RFQT).
}
\value{
\code{getPredict} returns a matrix where the row corresponds to the sample size (either the OOB data or the testing data, depending on the indicator used) and the column coresponds to the number of Q-trees. Each column represents the individual predicted effects with the present number of Q trees.
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

predict_matrix<-getPredict(RES,indicator=1 )
predict_matrix<-getPredict(RES,indicator=2 )
 
dim(predict_matrix)
View(predict_matrix)
}

