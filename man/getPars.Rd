\name{getPars}
\alias{getPars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get the present parameters information
}
\description{
\code{getPars} returns the information of the (hyper-)parameters relevant to the Q-tree and RFQT fitting.
}
\usage{
getPars(empty.argument)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{empty.argument}{
No arguments are needed.
}

}


\value{
\code{getPars()} returns the following information: the total number of candidate covariates, the training data size, the testing data siz, the stratification method, the size of pre-stratum (if applicable), the proportion of the covariates randomly considered in each node split, the maximal tree depth allowed, the way to calculate the instrument-exposure associations, the fixed instrument-exposure association level (if applicable), the minimal node size, the threshold value of the Q statistic.

If the variable is not defined, it will return 'Not defined' with a value in brackets representing the default value for this variable when fitting RFQT.

Note that only the variable defined in the global enviroment will be considered. Even if one defines these variables, the functions like \code{BootstrapTreeFitting} and \code{RFQTfit} will use the default value if (s)he does not modify the function arguments.

See the original paper for the details of each parameter information.
}

\author{
Haodong Tian
}

\examples{
library(MendelianRandomization)
library(tidyverse) 
library(data.table) 
library(parallel)

res<-getDat()
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set

method<-'DR'; SoP<-20; rate<-2/5

getPars()
}

