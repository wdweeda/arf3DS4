\name{makeExpDirs}
\Rdversion{1.1}
\alias{makeExpDirs}
\title{ Make Experiment directories. }
\description{
makeExpDirs makes an ARF experiment directory structure.
}
\usage{
makeExpDirs(path = getwd(), name = "default_experiment", subjectind = 1, 
 conditionind = 1, settings = new("settings"))
}
\arguments{
  \item{path}{Path in which to create the experiment directory.}
  \item{name}{Name of the experiment (=name of the directory).}
  \item{subjectind}{Vector of subject names for which to create directories.}
  \item{conditionind}{Vector of condition names for which to create directories.}
  \item{settings}{An optional "settings" object (see \code{\linkS4class{settings}}).}
}
\details{
makeExpDirs only makes the directory structure, data files must be copied to the directories manually (see \code{\linkS4class{data}}). After copying a call to \code{\link{loadExp}} in mode 'set', will set the experiment.Rda.
}
\value{
Returns an "experiment" object (\code{\linkS4class{experiment}}).
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{settings}}, \code{\link{loadExp}}, \code{\linkS4class{experiment}}
}
\keyword{ models }
