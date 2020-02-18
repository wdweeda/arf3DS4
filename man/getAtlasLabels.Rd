\name{getAtlasLabels}
\Rdversion{1.1}
\alias{getAtlasLabels}

\title{ Find anatomical atlas labels.}
\description{
getAtlasLabels finds anatomical atlas labels for a set of coordinates (either native ARF space or standard MNI space). It requires a working copy of FSL.
}
\usage{
getAtlasLabels(coordinates, registration, coortype = c("arf", "mni"), 
 atlas = c("both", "Talairach", "HarvardOxford"), ...)
}
\arguments{
  \item{coordinates}{Matrix with three columns (x,y,z) containing coordinates of a region. Multiple regions are in the rows of the matrix.}
  \item{registration}{A "registration" object belonging to the data the coordinates are from (see \code{\linkS4class{registration}}).}
  \item{coortype}{In what coordinate space are the coordinates (ARF native or MNI space).}
  \item{atlas}{For which atlases should the labels be returned? Currently only 'talairach' and 'harvardOxford' are available.}
  \item{\dots}{Arguments passed to the atlas lookup functions (\code{FSLDIR} indicates the FSL directory, \code{which} is '1mm' or '2mm')}
}
\details{
getAtlasLabels only works with a working copy of a full FSL installation (>=4.1).
}
\value{
Returns a list of length \code{nrow(coordinates)} with a \code{talairach} and/or a \code{harvard} slot containing the labels.
}
\references{
See \url{http://www.fmrib.ox.ac.uk/fsl/} on how to obtain FSL and licences. See \url{http://www.fmrib.ox.ac.uk/fsl/data/atlas-descriptions.html} on information on the atlases.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{registration}}, \code{\linkS4class{model}}, \code{\link{setRegFiles}}
}
\keyword{ models }
