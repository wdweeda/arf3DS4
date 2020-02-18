\name{checkSolution}
\alias{checkSolution}

\title{ Check the integrity of an ARF model}
\description{
Check the integrity of an ARF model for bounded parameters.}
\usage{
checkSolution(arfmodel, options = loadOptions(arfmodel), 
dat = readData(.model.avgdatfile(arfmodel)), thres = 6)
}

\arguments{
  \item{arfmodel}{A model-object (class \code{\linkS4class{model}})}
  \item{options}{An options-object (class \code{\linkS4class{options}})}
  \item{dat}{Average datafile.}
  \item{thres}{Number of digits for equality check.}
}
\value{
Returns the model object with added warnings if parameters are on a bound.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\keyword{ models }
