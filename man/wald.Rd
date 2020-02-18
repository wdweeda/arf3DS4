\name{wald}
\Rdversion{1.1}
\alias{wald}
\title{Calculate Wald statistics for an arf model.}
\description{Calculate Wald statistics for an arf model.}
\usage{
wald(arfmodel, waldobject = new("wald"), options = loadOptions(arfmodel))
}
\arguments{
  \item{arfmodel}{A valid model object.}
  \item{waldobject}{An object of class "wald". By default a new (empty) "wald" object is created. To test specific hypotheses a "wald" object with the \code{@consts} modified can be used.}
  \item{options}{An object of class "options".}
}
\details{
\code{wald} requires the (co)variance matrix of parameter estimates, it will give an error when the \code{arfmodel@varcov} slot does not exist (or has invalid values). 
}
\value{
Returns an object of class "model" with the appropriate slots.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}, \code{\link{varcov}}, \code{\link{processModel}}, \code{\linkS4class{options}}
}
\keyword{ models }
