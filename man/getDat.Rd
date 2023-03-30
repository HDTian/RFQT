\name{getDat}
\alias{getDat}
\title{
Creat a Toy Data
}
\description{
\code{getDat} creats a toy data based on a certain model
}
\usage{
getDat(N = 150000,
      Nt = 100000,
      Nc = 20,
      SoM = 0.5,
      ZXeffect = 0.5,
      XMeffect = 0.5,
      split = TRUE)
}
\arguments{
  \item{N}{
a integer indicates the total sample size (training data size + testing data size)
}
  \item{Nt}{
a integer indicates the size of the training data.
}
  \item{Nc}{
a integer indicates the number of candidate variables.
}
  \item{SoM}{
a value indicates the strength of modification (i.e. how strong will the covariate modify the treatment effect). Note even if \code{SoM=0}, there still exists weak modification.
}
  \item{ZXeffect}{
a value indicates the instrument-exposure effect.
}
  \item{XMeffect}{
a value indicates the effect of the exposure on the covariates.
}
  \item{split}{
logical. If \code{FALSE}, only one data set is returened. The default is \code{TRUE}, the data will be splitted into the training data and testing data.
}
}
\details{
The data-generating model is the Scenario B in the original paper, where the first 5 covariates are the effect modifiers. The covariates in even positions (2,4,...) are the common downstream variables of the exposure and te confounders, therefore causing collider bias (when conditioning on these variables).
}
\value{
\code{getDat()} returns a list with three kings of toy data set
\item{whole.data}{the complete whole dataset}
\item{traning.set}{the training set}
\item{testing.set}{the testing set}

the training set and the testing set are of the same form that can be regonized by all the functions. The first four columns are the individual IDs, the instrument, the exposure, and the outcome, respectively. The high-dimensional (dimensions = \code{Nc}) covariates are then follwing. The end column is \code{true_STE}, representing the individal controlled direct treatment effect (see more details in the original paper).
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Haodong Tian
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
res<-getDat()
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set
}

