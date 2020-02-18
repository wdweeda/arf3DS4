\name{loadFunc}
\Rdversion{1.1}
\alias{loadFunc}

\title{
Load functional data for a single run of a condition.
}
\description{
Load functional data (linking the run to the raw time-series) for a single run of a condition. This information is used by \code{\link{makeSingleTrialEvents}} to estimate single-trial data.
}
\usage{
loadFunc(subject, condition, run, experiment = NULL)
}

\arguments{
  \item{subject}{Subject name.}
  \item{condition}{Condition name.}
  \item{run}{Which run (can be numeric, or a directory name).}
  \item{experiment}{Optional experiment object.}
}
\value{
Returns an object of class \code{\linkS4class{functional}}.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{functional}}, \code{\linkS4class{experiment}}, \code{\link{setFuncTimings}}
}
\keyword{ models }