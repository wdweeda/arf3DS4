\name{fitRange}
\alias{fitRange}
\title{Fit a range of models.}
\description{fitRange performs model fits for a range of models.}
\usage{
fitRange(subject, condition, range, options=new('options'), 
 modelprefix='searchmodel', modeltype=c('gauss','simple'), 
 experiment=NULL)
}
\arguments{
  \item{subject}{Name of the subject.}
  \item{condition}{Name of the condition.}
  \item{range}{Vector containing the different numbers of regions for the spatial model.}
  \item{options}{An options object.}
  \item{modelprefix}{Prefix to append before the model name.}
  \item{modeltype}{Fit a full Gaussian or simple Gaussian model.}
  \item{experiment}{Optional experiment object (by default uses \code{.experiment}).}
}

\value{Returns an object of class \code{\linkS4class{sequence}} holding information on the valid models (minimum, BIC value, number of regions).}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}, \code{\linkS4class{options}}, \code{\linkS4class{experiment}}
}
\keyword{ models }