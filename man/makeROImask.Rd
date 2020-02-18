\name{makeROImask}
\alias{makeROImask}
\title{ Mask an fMRI image }
\description{
MakeROImask masks an fMRI image based on the non-zero elements of a maskfile (of class \code{\linkS4class{fmri.data}}).
}
\usage{
makeROImask(fmridata, maskdata)
}

\arguments{
  \item{fmridata}{The data to be masked (class \code{\linkS4class{fmri.data}}).}
  \item{maskdata}{The mask to be applied (voxels with non-zero-values are masked), can be of class \code{\linkS4class{fmri.data}} or a numeric vector.}
}
\value{Returns a masked object of class \code{\linkS4class{fmri.data}}.}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{fmri.data}}
}
\keyword{ models }