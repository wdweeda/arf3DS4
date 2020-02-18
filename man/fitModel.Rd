\name{fitModel}
\Rdversion{1.1}
\alias{fitModel}
\alias{fitModelOptim}
\alias{fitSimpleModelOptim}

\title{ Fit an ARF model.}
\description{
fitModel performs the fit-procedure of an ARF model to the data. It uses the values in the \code{\linkS4class{options}} object to modify the behavior of the fit-procedure.
By default \code{\link{optim}} (L-BFGS-B) is used for minimization. fitModel is a wrapper for \code{fitModelOptim} and \code{fitSimpleModelOptim}  
}
\usage{
fitModel(arfmodel, options = loadOptions(arfmodel), 
 dat = readData(.model.avgdatfile(arfmodel)), 
 weights = readData(.model.avgWfile(arfmodel)), 
 printlevel = 0, try.silen = T)
fitModelOptim(arfmodel, options = loadOptions(arfmodel), 
 dat = readData(.model.avgdatfile(arfmodel)), 
 weights = readData(.model.avgWfile(arfmodel)), 
 printlevel = 0, try.silen = T)
fitSimpleModelOptim(arfmodel, options = loadOptions(arfmodel), 
 dat = readData(.model.avgdatfile(arfmodel)), 
 weights = readData(.model.avgWfile(arfmodel)), 
 printlevel = 0, try.silen = T)

}
\arguments{
  \item{arfmodel}{An arf model object (class \code{\linkS4class{model}})}
  \item{options}{An "options" object (class \code{\linkS4class{options}}).The \code{options@start.method} slot is used to determine which method to use for starting values. The \code{arfmodel@modeltype} slot is used to determine whether to fit a full Gaussian model (default) or a simple Gaussian model (with fixed width parameters).}  
  \item{dat}{An object of class "fmri.data" (\code{\linkS4class{fmri.data}} containing the averaged beta-values.)}
  \item{weights}{An object of class "fmri.data" (\code{\linkS4class{fmri.data}} containing the averaged weights of the beta-values.)}
  \item{printlevel}{Number indicating the level of output of the minimization procedure.}
  \item{try.silen}{Should errors be printed when they occur or only passed to the \code{warnings} slot.}
}
\details{
fitModel calls one of the following functions based on the "options" object: \code{fitModelOptim}, \code{fitSimpleModelOptim}. These functions can also be called directly. Models can be created by calling \code{\link{newModel}}.
}
\value{
Returns an object of class "model" with the appropriate slots.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}, \code{\linkS4class{options}}, \code{\linkS4class{fmri.data}}, \code{\link{newModel}}
}
\keyword{ models }
