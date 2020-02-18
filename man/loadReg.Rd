\name{loadReg}
\Rdversion{1.1}
\alias{loadReg}
\title{Load a registration object.}
\description{
Load a registration object for a single run of a condition.
}
\usage{
loadReg(subject, condition, run, experiment = NULL)
}
\arguments{
  \item{subject}{Subject.}
  \item{condition}{Condition.}
  \item{run}{Which run (can be numeric, or a directory name)}
  \item{experiment}{Optional experiment object.}
}
\value{
Returns an object of class \code{\linkS4class{registration}}.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{registration}}, \code{\linkS4class{experiment}}
}
\keyword{ models }
