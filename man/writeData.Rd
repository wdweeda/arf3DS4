\name{writeData}
\Rdversion{1.1}
\alias{writeData}
\title{Write NIFTI/ANALYZE files.}
\description{Write a data-vector to a NIFTI/ANALYZE file.}
\usage{
writeData(headinf, datavec)
}
\arguments{
  \item{headinf}{An object of class "fmri.data".} 
  \item{datavec}{An optional data-vector to be written to the nifti/analyze file. If no \code{datavec} is given the \code{@datavector} slot of the fmri.data-object is used.}
}
\details{
All information on filetype \code{@filetype}, path \code{@fullpath}, and filename \code{@filename} are within the fmri.data-object (also see \code{\linkS4class{nifti.fileinfo}}). 
Changing these values will modify where and how the data are saved. 
}
\value{
Returns \code{TRUE} when successful.
}
\references{ See the Nifti homepage (\url{http://nifti.nimh.nih.gov/}) for more information on the nifti filetype.}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{fmri.data}}, \code{\linkS4class{nifti.header}}, \code{\linkS4class{nifti.fileinfo}}, \code{\link{readData}}
}
\keyword{ models }
\keyword{ methods }

