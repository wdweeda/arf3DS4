\name{setRegFiles}
\Rdversion{1.1}
\alias{setRegFiles}
\title{Set registration files.}
\description{
Link the files necessary for registration to the appropriate slots in the registration object. Registration requires a working copy of FSL. Files must be copied from the FSL directory to the \bold{/registration} directory of the ARF experiment.
By default the standard filenames used by FSL are used. 
}
\usage{
setRegFiles(registration, examp2stand = "example_func2standard.mat", 
 examp2high = "example_func2highres.mat", high2stand = "highres2standard.mat", 
 example_func = "example_func.nii.gz", highres = "highres.nii.gz", 
 standard = "standard.nii.gz")
}
\arguments{
  \item{registration}{Object of class \code{\linkS4class{registration}} to link to.}
  \item{examp2stand}{ARF to standard space affine transformation matrix file. }
  \item{examp2high}{ARF to high resolution affine transformation matrix file. }
  \item{high2stand}{High resolution to standard space affine transformation matrix file.}
  \item{example_func}{Nifti file in ARF native space. }
  \item{highres}{T1 weighted structural image.}
  \item{standard}{Standard space image (MNI152).}
}
\value{
Returns a registration object with the appropriate slots filled.
}
\references{
See \url{http://www.fmrib.ox.ac.uk/fsl/} on how to obtain FSL and licences.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{registration}}, \code{\link{createRegs}}, \code{\link{setRegParams}}
}
\keyword{ methods }
\keyword{ models }

