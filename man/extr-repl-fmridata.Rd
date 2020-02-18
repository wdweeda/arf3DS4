\name{extr-repl-fmridata}
\alias{[,fmri.data-method}
\alias{[<-,fmri.data-method}
\title{Extract/replace elements of an fmri.data object}
\description{Extract or replace elements of an fmri.data object}
\details{fMRI data can be accessed using standard R array indexing. The dimensions of the array are determined by the dimensions of the \code{\linkS4class{fmri.data}} object. Both 3D and 4D (time-series) data can be used.}
\value{Returns an object of class "fmri.data" (in the case of replacement), with the appropriate values replaced) Or returns the appropriate part of the array.}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{fmri.data}}, \code{\link{fmri2array}}
}
\keyword{ models }