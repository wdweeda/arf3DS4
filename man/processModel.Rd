\name{processModel}
\Rdversion{1.1}
\alias{processModel}
\title{Process an ARF model.}
\description{Run an entire processing pipeline for a given ARF model (minimization and statistics).}
\usage{
processModel(arfmodel, options = loadOptions(arfmodel), 
 dat = readData(.model.avgdatfile(arfmodel)), 
 weights = readData(.model.avgWfile(arfmodel)), 
 pr = T, printlevel = 0, try.silen = T)
}
\arguments{
  \item{arfmodel}{An arf model object (class \code{\linkS4class{model}}).}
  \item{options}{An "options" object (class \code{\linkS4class{options}}).The \code{options@start.method} slot is used to determine which method to use for starting values. 
  The \code{arfmodel@modeltype} slot is used to determine whether to fit a full Gaussian model (default) or a simple Gaussian model.}  
  \item{dat}{An object of class "fmri.data" (\code{\linkS4class{fmri.data}} containing the averaged beta-values.)}
  \item{weights}{An object of class "fmri.data" (\code{\linkS4class{fmri.data}} containing the averaged weights of the beta-values.)}
  \item{pr}{Print output to the screen, or perform silently.}
  \item{printlevel}{Number indicating the level of output of the minimization procedure.}
  \item{try.silen}{Should errors be printed when they occur or only passed to the \code{@warnings} slot.}
}
\details{
processModel calls \code{\link{fitModel}}, \code{\link{varcov}} and \code{\link{wald}}. The model is automatically saved.
}
\value{
Returns an object of class "model" with the appropriate slots.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}, \code{\linkS4class{options}}, \code{\linkS4class{fmri.data}}, \code{\link{newModel}}, \code{\link{fitModel}}, \code{\link{varcov}}, \code{\link{wald}}
}
\keyword{ models }
