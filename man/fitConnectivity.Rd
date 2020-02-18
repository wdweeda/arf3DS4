\name{fitConnectivity}
\alias{fitConnectivity}

\title{
Estimate connectivity of an ARF model
}
\description{
fitConnectivity estimates the functional connectivity between activated regions in an ARF model based on trial-to-trial amplitude variation. 
It uses the single-trial data (extracted from the raw time-series using \code{\link{makeSingleTrialEvents}}), to estimate the trial-by-trial amplitudes of
each region in the ARF model.
}
\usage{
fitConnectivity(arfmodel, funcfilename="single_events.nii.gz")
}

\arguments{
  \item{arfmodel}{An ARF model object (class \code{\linkS4class{model}})}
  \item{funcfilename}{Filename of the volume containing the single-trial data.}
}
\details{
Each volume of the file specified in funcfilename should contain estimates of a single-trial (see \code{\link{makeSingleTrialEvents}}).
}
\value{
Returns an object of class \code{\linkS4class{arfcorrelation}} including the region by time matrix of amplitude time-series and the estimated correlations.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}, \code{\link{makeSingleTrialEvents}}, \code{\linkS4class{arfcorrelation}}, \code{\link{roiConnectivity}}
}
\keyword{ models }
