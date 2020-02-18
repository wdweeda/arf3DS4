\name{setIsoContour}
\alias{setIsoContour}

\title{Create isocontours}
\description{Create isocontours for each region in the ARF model.}
\usage{
setIsoContour(arfmodel, conf.int = 95)
}
\arguments{
  \item{arfmodel}{An ARF model object (class \code{\linkS4class{model}}).}
  \item{conf.int}{Confidence Interval (in percentage) for the isocontour.}
}
\value{An object of class \code{\linkS4class{fmri.data}} with a volume for each region in the spatial model. Each volume indicates the voxels within the isocontour for that specific region.}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}, \code{\linkS4class{fmri.data}}, \code{\link{roiConnectivity}}
}
\keyword{ models }
