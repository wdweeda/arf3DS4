\name{getModelAtlas}
\Rdversion{1.1}
\alias{getModelAtlas}
\title{ Get Atlas Labels for an ARF model.}
\description{
getModelAtlas fetches atlas labels for the center coordinates of all regions in a model.
}
\usage{
getModelAtlas(arfmodel, regrun = 1, saveastext = F)
}
\arguments{
  \item{arfmodel}{Object of class "model" (see \code{\linkS4class{model}}).}
  \item{regrun}{Which run is used for the coordinate conversion?}
  \item{saveastext}{Save a text file in the working directory with atlas information?}
}
\details{
getModelAtlas only works with a working copy of a full FSL installation (>=4.1). It calls \code{\link{getAtlasLabels}}.
}
\value{
Returns a list of length \code{nrow(coordinates)} with a \code{talairach} and/or a \code{harvard} slot containing the labels.
}
\references{
See \url{http://www.fmrib.ox.ac.uk/fsl/} on how to obtain FSL and licences. See \url{http://www.fmrib.ox.ac.uk/fsl/data/atlas-descriptions.html} on information on the atlases.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}, \code{\link{getAtlasLabels}}
}
\keyword{ models }
