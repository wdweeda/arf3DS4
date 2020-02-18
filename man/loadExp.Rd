\name{loadExp}
\Rdversion{1.1}
\alias{loadExp}
\title{ Load an arf3DS4 experiment }
\description{
loadExp loads an arf3DS4 experiment.
}
\usage{
loadExp(path = getwd(), method = c("fast", "set", "rda"))
}
\arguments{
  \item{path}{Path to the experiment directory.}
  \item{method}{Which method is used for loading.}
}
\details{
If \code{method='fast'} (default) experiment.Rda is loaded without checks directly.
If \code{method='set'} the experiment is entirely build up by searching the directories of the experiment for subjects/conditions and models. Use this method if the experiment has changed (added subjects, moved to another directory etc.).
If \code{method='rda'} the experiment is checked using the data in the experiment.Rda (except for the root-path). Use this method if the experiment has been moved (but the subjects/conditions have not changed).
loadExp creates an object named \code{.experiment} in the \code{.arfInternal} environment, this is to facilitate calls to functions which have an 'experiment' object as input. All ARF functions get the \code{.experiment} object from \code{.arfInternal}.
}
\value{
Returns an object of class \code{\linkS4class{experiment}} invisibly.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{experiment}}
}
\keyword{ models }
