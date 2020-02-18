\name{loadData}
\Rdversion{1.1}
\alias{loadData}

\title{Load an arf3DS4 "data" object.}
\description{
loadData loads a "data" object for a given subject and condition.
}
\usage{
loadData(subject, condition, experiment = NULL)
}
\arguments{
  \item{subject}{Name of the subject.}
  \item{condition}{Name of the condition.}
  \item{experiment}{Optional object of class "experiment" (see \code{\linkS4class{experiment}})}
}
\details{
loadData uses the loaded experiment as a default (see \code{\link{loadExp}}) on how to load experiments.
}
\value{
Returns an object of class "data" (see \code{\linkS4class{data}}).
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{data}}, \code{\linkS4class{experiment}}, \code{\link{loadExp}}
}
\keyword{ models }

