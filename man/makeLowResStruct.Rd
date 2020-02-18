\name{makeLowResStruct}
\Rdversion{1.1}
\alias{makeLowResStruct}
\alias{makeLowResStructAvg}
\title{Make low resolution structural image from high resolution T1 image.}
\description{Make low resolution structural image from high resolution T1 image. The function requires that registration files (\code{\link{setRegFiles}}) and parameters (\code{\link{setRegParams}}) are set. It also requires a working copy of FSL.}
\usage{
makeLowResStruct(arfdata, experiment = NULL)
makeLowResStructAvg(arfmodel, experiment = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{arfdata}{An arf "data" object (see \code{\linkS4class{data}}.)}
  \item{arfmodel}{An arf "model" object (see \code{\linkS4class{model}}.)}
  \item{experiment}{An arf "experiment" object (see \code{\linkS4class{experiment}}.)}
}
\details{
makeLowResStructAvg makes use of the images created by makeLowResStruct so this must be called first. makeLowResStruct has to be run only once for each condition. 
}
\value{
None.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\references{
See \url{http://www.fmrib.ox.ac.uk/fsl/} on how to obtain FSL and licences.
}

\seealso{
	\code{\linkS4class{model}}, \code{\linkS4class{experiment}}, \code{\linkS4class{data}}, \code{\link{setRegFiles}}, \code{\link{setRegParams}}
}
\keyword{ models }
