\name{getDat}
\alias{getDat}
\title{
Creat a Simulated Data
}
\description{
\code{getDat} creats a toy data based on a certain model
}
\usage{
getDat(N = 150000,
      Nt = 100000,
      Nc = 20,
      scenario='A',
      SoM = 0.5,
      ZXeffect = 0.5,
      Random = TRUE,
      label = TRUE,
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
  \item{Scenario}{
a character indicates the model scenarios data will be simulated from. \code{'A'} refers to scenario where no colliders. \code{'B'} refers to the scenario that half covariates will be colliders. \code{'C'} refers to more complicated collider scenarios.
}
  \item{SoM}{
a value indicates the strength of modification (i.e. how strong will the covariate modify the treatment effect). Note even if \code{SoM=0}, there still exists weak modification.
}
  \item{ZXeffect}{
a value indicates the instrument-exposure effect.
}
  \item{Random}{
logical. Indicates whether to use the random value for strenght modification. Default value is \code{TRUE}, where weak modification always exists even if \code{SoM=0} due to randomness.  
}
 \item{label}{
logical. Indicates whether to added the final column of the simulated data as the true heterogenerous effect (i.e. label). Default value is \code{TRUE}.
}
  \item{split}{
logical. If \code{FALSE}, only one data set is returened. The default is \code{TRUE}, the data will be splitted into the training data and testing data.
}
}
\details{
The data-generating model is the same as the Scenarios in the original paper, where the first 5 covariates are the true effect modifiers. The covariates in even positions (2,4,...) are the common downstream variables of the exposure and te confounders, therefore causing collider bias (when conditioning on these variables).
}
\value{
\code{getDat()} returns a list with three kings of toy data set
\item{whole.data}{the complete whole dataset}
\item{traning.set}{the training set}
\item{testing.set}{the testing set}
\item{SoM}{Strength of Modification}
\item{modifier_vec}{the modification strength vector for each covariate}
\item{Scenario}{scenario information}

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
set.seed(60)
res<-getDat() #simulated data #the deflaut setting: scenario='A' and SoM=0.5
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set
}

