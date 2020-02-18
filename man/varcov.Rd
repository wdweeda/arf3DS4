\name{varcov}
\Rdversion{1.1}
\alias{varcov}
\title{Calculate (co)variance matrix of an ARF model.}
\description{Calculate (co)variance matrix of the parameters in an ARF model using Sandwich estimation. The method of Sandwich estimation can be modified via the \code{\linkS4class{options}} object of the ARF model.}
\usage{
varcov(arfmodel)
}
\arguments{
  \item{arfmodel}{A valid ARF model object.}
}
\details{
\code{varcov} needs first-order derivatives and residual matrices, if they do not exist, they are created automatically.
}
\value{
Returns an object of class "model" with the appropriate slots.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}, \code{\link{processModel}}, \code{\link{makeDerivs}}, \code{\linkS4class{options}} 
}
\keyword{ models }