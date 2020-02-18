\name{minBIC}
\alias{minBIC}
\title{ Calculate BIC values for multiple ARF models. }
\description{Calculate and show the BIC, minimum and number of regions of all models in a condition.}
\usage{
minBIC(subject, condition)
}
\arguments{
  \item{subject}{Name of the subject.}
  \item{condition}{Name of the condition.}
}

\value{Returns an object of class \code{\linkS4class{sequence}} containing fit-information of all models.}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{fmri.data}}, \code{\linkS4class{sequence}}
}
\keyword{ models }
