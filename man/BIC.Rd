\name{BIC}
\Rdversion{1.1}
\alias{BIC}
\title{Bayesian Information Criterion}
\description{
BIC calculates the Bayesian Information Criterion (model fit) for a valid ARF model. BIC requires a valid (converged) ARF model (\code{\linkS4class{model}}). }
\usage{
BIC(arfmodel, options = loadOptions(arfmodel))
}

\arguments{
  \item{arfmodel}{A fitted (and valid) ARF model object of class \code{\linkS4class{model}}. }
  \item{options}{An "options" object. The options object belonging to the ARF model is loaded by default.}
}

\details{
BIC calculates the BIC value using the number of brain voxels (Set by the \code{model@mask} slot). Without constants the BIC equals:

	\code{BIC = ln(S(theta)) + p * ln(n)} 

where \code{p} is the number of parameters in the spatial model, \code{n} is the number of brain voxels, and \code{S(theta)} denotes the minimum of the objective function.


}

\value{
Returns an object of class "model" with the slot \code{model@fit[1]} set to the BIC value. 
}

\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\link{RMSEA}}, \code{\linkS4class{model}}, \code{\linkS4class{options}}
}
\keyword{ models }
