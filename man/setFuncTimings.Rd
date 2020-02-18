\name{setFuncTimings}
\Rdversion{1.1}
\alias{setFuncTimings}
\title{Set the timings of the stimuli and link to raw time-series.}
\description{
Set the timings of the stimuli (i.e. single trials) and link these timings to raw time-series. For each run of each condition the timinigs of when trials were presented
must be given for ARF to determine the single-trial data.
}
\usage{
setFuncTimings(subject, condition, run, timings, func_data = NULL, 
 experiment = NULL)
}
\arguments{
  \item{subject}{Name of the subject.}
  \item{condition}{Name of the condition.}
  \item{run}{Which run, can be numeric or a filename.}
  \item{timings}{A vector with timings (in seconds), possibly with an attribute \code{stimlen} (vector of stimulus lengths).}
  \item{func_data}{Name of the raw time-series datafile (located in the \bold{/subject/funcs} directory).}
  \item{experiment}{An optional experiment object.}
}
\details{
An object of class \code{\linkS4class{functional}} is also saved in the appropriate directory.
}
\value{
Returns an object of class \code{\linkS4class{functional}}.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{functional}}, \code{\linkS4class{experiment}}, \code{\link{makeSingleTrialEvents}}
}
\keyword{ methods }
\keyword{ models }
