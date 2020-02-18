\name{RMSEA}
\Rdversion{1.1}
\alias{RMSEA}
\title{Calculate Root Mean Square Error.}
\description{
Calculate RMSEA fit value.
}
\usage{
RMSEA(arfmodel, options = loadOptions(arfmodel))
}
\arguments{
  \item{arfmodel}{A fitted (and valid) arf-model object of class "model". }
  \item{options}{An "options" object. The options object belonging to the ARF model is loaded by default.}
}

\details{
RMSEA calculates the RMSEA value using the number of brain voxels (Set by the \code{model@mask} slot).
}

\value{
Returns an object of class "model" with the slot \code{model@fit[2]} set to the RMSEA value. 
}

\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\link{BIC}}, \code{\linkS4class{model}}, \code{\linkS4class{options}}
}
\keyword{ models }

