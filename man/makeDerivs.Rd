\name{makeDerivs}
\Rdversion{1.1}
\alias{makeDerivs}
\title{ Calculate analytical first-order derivatives of parameters in an ARF model.}
\description{
makeDerivs calculates analytical first-order derivatives of the parameters in an ARF model. It writes the derivative matrix to a binary file.
}
\usage{
makeDerivs(arfmodel,method=c('viaR','direct'))
}
\arguments{
  \item{arfmodel}{An arf model object (see \code{\linkS4class{model}}).}
  \item{method}{Specification of derivative saving method. \code{'viaR'} is the default.}
}
\value{
Returns \code{TRUE} when successful.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}
}
\keyword{ models }
