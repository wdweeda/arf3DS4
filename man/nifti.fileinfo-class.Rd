\name{nifti.fileinfo-class}
\Rdversion{1.1}
\docType{class}
\alias{nifti.fileinfo-class}
\alias{.nifti.fileinfo.endian,nifti.fileinfo-method}
\alias{.nifti.fileinfo.endian<-,nifti.fileinfo-method}
\alias{.nifti.fileinfo.extension,nifti.fileinfo-method}
\alias{.nifti.fileinfo.extension<-,nifti.fileinfo-method}
\alias{.nifti.fileinfo.filename,nifti.fileinfo-method}
\alias{.nifti.fileinfo.filename<-,nifti.fileinfo-method}
\alias{.nifti.fileinfo.filetype,nifti.fileinfo-method}
\alias{.nifti.fileinfo.filetype<-,nifti.fileinfo-method}
\alias{.nifti.fileinfo.fullpath,nifti.fileinfo-method}
\alias{.nifti.fileinfo.fullpath<-,nifti.fileinfo-method}
\alias{.nifti.fileinfo.gzipped,nifti.fileinfo-method}
\alias{.nifti.fileinfo.gzipped<-,nifti.fileinfo-method}
\alias{.nifti.fileinfo.version,nifti.fileinfo-method}
\alias{.nifti.fileinfo.version<-,nifti.fileinfo-method}

\title{arf3DS4 "nifti.fileinfo" class }
\description{The arf3DS4 "nifti.fileinfo" class contains the basic file information needed to read/write an ANALYZE/NIFTI header.}
\section{Objects}{
Objects can be created by calls of the form \code{new("nifti.fileinfo", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{fullpath}:}{Full path of fMRI datafile.}
    \item{\code{filename}:}{Filename of the fMRI datafile. }
    \item{\code{filetype}:}{Type of the fMRI datafile. }
    \item{\code{extension}:}{File extension. }
    \item{\code{gzipped}:}{Is the file gzipped? }
    \item{\code{endian}:}{Endianness of file. }
    \item{\code{version}:}{Object of class "version" (see \code{\linkS4class{version}})}
  }
}
\section{Methods}{
 Standard arf3DS4 accessor and replacement functions can be used to access and replace slots. 
  \describe{
 	\item{access}{\code{.classname.slotname(object)}}
    \item{replace}{\code{.classname.slotname(object)<-}}
 	 }
 In this case classname is "nifti.fileinfo" and slotname is for example "extension". So to access the "extension" slot of an r-object "myniftifileinfo"  (of class "nifti.fileinfo") type \code{.nifti.fileinfo.extension(myfniftifileinfo)}.
}
\references{ See the Nifti homepage (\url{http://nifti.nimh.nih.gov/}) for more information on the nifti filetype.}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{nifti.header}}, \code{\linkS4class{fmri.data}}
}
\examples{
showClass("nifti.fileinfo")
}
\keyword{classes}
