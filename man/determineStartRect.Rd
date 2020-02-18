\name{determineStartRect}
\Rdversion{1.1}
\alias{determineStartRect}
\alias{determineStartRectSimple}

\title{Calculate starting values for an ARF model}
\description{
determineStartRect calculates starting values for an ARF model by searching for maxima and minima in the data, and estimating the width of these maxima/minima.
}
\usage{
determineStartRect(arfmodel, options = loadOptions(arfmodel))
determineStartRectSimple(arfmodel, options = loadOptions(arfmodel))
}
\arguments{
  \item{arfmodel}{An ARF model object (class \code{\linkS4class{model}})}
  \item{options}{An "options" object (class \code{\linkS4class{options}}), by default the options.Rda from the current model is loaded. The value of \code{options@start.maxfac} determines the calculation of the spatial extent values, higher settings will lead to greater extent values.}
}
\details{
determineStartRect determines the maxima and minima of the data and at each maximum/minimum searches in each dimension (x,y,z) the point where the maximum/minimum value is half its original value. 
This creates a box around the maximum/minimum, which is set to zero, to avoid searching the same area twice. The widths of this box are taken as starting values for the width-parameters.
}
\value{
Returns the ARF "model" object, with the starting values in the \code{@startval} slot.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}, \code{\linkS4class{options}}
}
\keyword{ models }
