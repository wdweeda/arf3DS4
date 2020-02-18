\name{setRegParams}
\Rdversion{1.1}
\alias{setRegParams}
\title{Set the registration parameters.}
\description{Set the registration matrices to the appropriate values given the registration files. Registration requires a working copy of FSL.}
\usage{
setRegParams(registration)
}
\arguments{
  \item{registration}{Object of class \code{\linkS4class{registration}} with appropriate filenames.}
}
\details{
First run \code{\link{setRegFiles}} to set the correct filenames.
}
\value{
Returns a registration object with the appropriate slots filled.
}
\references{
See \url{http://www.fmrib.ox.ac.uk/fsl/} on how to obtain FSL and licences.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{registration}}, \code{\link{setRegFiles}}, \code{\link{createRegs}}
}
\keyword{ methods }
\keyword{ models }
