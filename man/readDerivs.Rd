\name{readDerivs}
\Rdversion{1.1}
\alias{readDerivs}
\title{Read (binary) derivative file}
\description{Read in derivatives of an ARF model.}
\usage{
readDerivs(arfmodel)
}
\arguments{
  \item{arfmodel}{A valid ARF model object.}
}
\value{
Returns an object of class "matrix" containing the derivatives (voxels x parameters).
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}, 
}
\keyword{ models }