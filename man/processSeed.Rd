\name{processSeed}
\alias{processSeed}
\title{Process a model using a seed number of regions}
\description{
processSeed runs a sequence of models with different starting values to obtain a stable solution. It starts with as simple model with rectangular starts, 
subsequently it runs a full model given the estimates of the simple model. It then prunes the full model, discarding the regions with bad gradients and/or 
bounded parameters.
}
\usage{
processSeed(modelname='defaultmodel', seedreg, subject='', 
 condition='', startmethod=c('default','simple'), grad=NULL, 
 bound=NULL, pval=NULL, options=new('options'), pr=T, 
 printlevel=0, try.silen=T, overwrite=T, experiment=NULL)
}

\arguments{
  \item{modelname}{Name to give the pruned models in the sequence.}
  \item{seedreg}{Number of regions to fit at sequence start (seed number).}
  \item{subject}{Name of the subject.}
  \item{condition}{Name of the condition.}
  \item{startmethod}{Use simple model as start or the default (rectangular method).}
  \item{grad}{Gradient check limits, absolute higher values are pruned.}
  \item{bound}{Threshold to check if a parameter is on a bound. When NULL, check is suppressed.}
  \item{pval}{Prune regions with non-significant Wald statistics higher than \code{pval}.}
  \item{options}{Options object.}
  \item{pr}{Output to screen?}
  \item{printlevel}{Number indicating the level of output of the minimization procedure.}
  \item{try.silen}{Should errors be printed when they occur or only passed to the \code{warnings} slot.}
  \item{overwrite}{Overwrite the model if a model with the same name already exists?}
  \item{experiment}{Optional experiment object.}
}

\value{
Returns the best model (class \code{\linkS4class{model}})
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}, \code{\link{newModel}}, \code{\link{fitModel}}, \code{\link{pruneModel}}, \code{\linkS4class{options}}, \code{\linkS4class{experiment}}
}
\keyword{ models }

