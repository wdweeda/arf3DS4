\name{euclidDist}
\Rdversion{1.1}
\alias{euclidDist}

\title{ Calculate Euclidian distances of all regions in a model. }
\description{
euclidDist calculates Euclidian distances between the center locations of each region in an ARF model. This can be used to assess if regions are close together.
}
\usage{
euclidDist(arfmodel, thres = 5, quiet = T)
}
\arguments{
  \item{arfmodel}{An arf model object (class \code{\linkS4class{model}})}
  \item{thres}{Threshold at which to output region information. If regions are closer than \code{thres} information is printed if \code{quiet=FALSE}.}
  \item{quiet}{Should region information be printed?}
}
\details{
euclidDist only uses center location, no information on the extent of a region is used.
}
\value{
Returns a lower triangular matrix with Euclidian distances of all regions in the model.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}
}
\keyword{ models }
