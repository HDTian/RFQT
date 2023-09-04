\name{GetNindex}
\alias{GetNindex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get the Q Tree End Node Index for Choosen Samples
}
\description{
\code{GetNindex} helps to match the Q-tree end node for any inputed samples, and therefore enables to make further analysis like effect prediction for these inputed samples.
}
\usage{
GetNindex(M, 
          rdat, 
          S = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A matrix or data frame contains the candidate covariate information.
}
  \item{rdat}{
A data frame. It is the Q-tree informaton result and obtained by \code{rdat<-GetTree()}.
}
  \item{S}{A positive integer indicates the max depth of each sample will explore along the tree. This number should be consistent with the max depth of the fitted Q tree (also refected in \code{rdat}). The default valus if \code{NA}, which means the code will automatically use the correct \code{S} value.
%%     ~~Describe \code{S} here~~
}
}
\details{
The inputed data set \code{M} must contain the same candidate covariate information (with the same order) as the training data for the Q-tree. In addition, the first four columns (individual IDs, the instrument, the exposure, the outcome) need to be removed, so that the fist column starts from the first candidate covariate.
}
\value{
\code{GetNindex} returns a vector with the same length of the row number of the inputed data \code{M}. Each element represents the end node (i.e. leaf) index according to the reference Q-tree (the tree information was stored within \code{rdat}).
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
res<-getDat() #simulated data #the deflaut setting: scenario='A' and SoM=0.5
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set
 
rdat<-GetTree(odat)
 
 
vdat_Nindex<-GetNindex(vdat[,5:24] ,rdat )#the M colnumber and order should be same as the training data
 
 
#may use another independent data (estimation data) to calculate the endnode-specific IV/MR estimates
 
#for example: let odat halves into trdata (tree data) and estdata (estimation data); also have rdat<-GetTree(trdata)
 
trdat<-odat[1: (nrow(odat)/2) ,]#tree data
estdat<-odat[(nrow(odat)/2+1):nrow(odat) ,] #estimaiton data
 
rdat<-GetTree(trdat)#fitted Q tree
 
estdat$Nindex<-GetNindex(estdat[,5:24] ,rdat ) #rdat contains the tree information (i.e. decision rule)
}
