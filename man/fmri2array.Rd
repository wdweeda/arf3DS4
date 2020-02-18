\name{fmri2array}
\alias{fmri2array}
\title{Convert fMRI data to an array}
\description{Convert the data vector of an fmri.data object to an array.}
\usage{
fmri2array(fmridat)
}
\arguments{
  \item{fmridat}{An object of class \code{\linkS4class{fmri.data}}.}
}
\value{Returns an array containing the fMRI data (\code{@datavec}-slot) with dimensions set by the fmri.data object.}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{fmri.data}}
}
\keyword{ models }