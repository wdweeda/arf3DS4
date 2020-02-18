\name{roiConnectivity}
\alias{roiConnectivity}

\title{calculate Connectivity with spatial averaging}
\description{Calculates functional connectivity based on the spatial average of the voxels within the isoContours of the estimated model.}
\usage{
roiConnectivity(arfmodel,roidata=setIsoContour(arfmodel,95),
 funcfilename='single_events.nii.gz',type=c('avg','ev'), 
 evmodel=c('spatial','spatiotemporal','eigenvariate'))
}

\arguments{
  \item{arfmodel}{An ARF model object (class \code{\linkS4class{model}}).}
  \item{roidata}{\code{\linkS4class{fmri.data}} object containing ROIs for each region (must contain a volume for each region), can be obtained via \code{\link{setIsoContour}}.}
  \item{funcfilename}{Filename of the volume to extract the single-trial data.}
  \item{type}{Average using no weighing (\code{'avg'}) or by weighing with the first spatial eigenvector (\code{'ev'}).}
  \item{evmodel}{Which data is used to calculate the first spatial eigenvector: \code{'spatial'}, \code{'spatiotemporal'}, or \code{'eigenvariate'}. } 
}

\value{
Returns an object of class \code{\linkS4class{arfcorrelation}} including the region by time matrix of amplitude time-series and the estimated correlations.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}, \code{\linkS4class{arfcorrelation}}, \code{\link{fitConnectivity}}, \code{\link{setIsoContour}}
}
\keyword{ models }