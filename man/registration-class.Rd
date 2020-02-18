\name{registration-class}
\Rdversion{1.1}
\docType{class}
\alias{registration-class}
\alias{.registration.Aex2hi,registration-method}
\alias{.registration.Aex2hi<-,registration-method}
\alias{.registration.Ahi2st,registration-method}
\alias{.registration.Ahi2st<-,registration-method}
\alias{.registration.Dex,registration-method}
\alias{.registration.Dex<-,registration-method}
\alias{.registration.Dhi,registration-method}
\alias{.registration.Dhi<-,registration-method}
\alias{.registration.Dst,registration-method}
\alias{.registration.Dst<-,registration-method}
\alias{.registration.examp2high,registration-method}
\alias{.registration.examp2high<-,registration-method}
\alias{.registration.examp2stand,registration-method}
\alias{.registration.examp2stand<-,registration-method}
\alias{.registration.example,registration-method}
\alias{.registration.example<-,registration-method}
\alias{.registration.filename,registration-method}
\alias{.registration.filename<-,registration-method}
\alias{.registration.fullpath,registration-method}
\alias{.registration.fullpath<-,registration-method}
\alias{.registration.high2stand,registration-method}
\alias{.registration.high2stand<-,registration-method}
\alias{.registration.highres,registration-method}
\alias{.registration.highres<-,registration-method}
\alias{.registration.linkedfile,registration-method}
\alias{.registration.linkedfile<-,registration-method}
\alias{.registration.OXst,registration-method}
\alias{.registration.OXst<-,registration-method}
\alias{.registration.standard,registration-method}
\alias{.registration.standard<-,registration-method}
\alias{.registration.SXhi,registration-method}
\alias{.registration.SXhi<-,registration-method}
\alias{.registration.version,registration-method}
\alias{.registration.version<-,registration-method}

\title{arf3DS4 "registration" class  }
\description{The arf3DS4 registration class contains information on the conversion matrices for brain registration purposes. 
It is used to convert ARF models to standard space (MNI152), and to structural (T1) images. It is also used to look up the location of activated regions in standard atlases (Harvard-Oxford and Talairach). 
The "registration" class is optional, for it to work properly a functioning installation of FSL is required.}
\section{Objects}{
Objects can be created by calls of the form \code{new("registration", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{fullpath}:}{Full path of registration directory.}
    \item{\code{filename}:}{Registration filename. }
    \item{\code{linkedfile}:}{Path to linked beta-file.}
    \item{\code{examp2high}:}{ARF to high resolution affine transformation. }
    \item{\code{high2stand}:}{High resolution to standard space affine transformation.}
    \item{\code{examp2stand}:}{ARF to standard space affine transformation. }
    \item{\code{example}:}{Nifti file in arf space. }
    \item{\code{highres}:}{T1 weighted structural image.}
    \item{\code{standard}:}{Standard space image (MNI152).}
    \item{\code{Dex}:}{ARF voxel to ARF mm. }
    \item{\code{Dhi}:}{highres mm to highres voxels. }
    \item{\code{Dst}:}{standard mm to standard voxels. }
    \item{\code{SXhi}:}{x-axis flip.}
    \item{\code{Aex2hi}:}{ARF mm to highres mm. }
    \item{\code{Ahi2st}:}{highres mm to standard mm. }
    \item{\code{OXst}:}{standard origin offset. }
    \item{\code{version}:}{Object of class "version" (see \code{\linkS4class{version}})}
  }
}
\section{Methods}{
 Standard arf3DS4 accessor and replacement functions can be used to access and replace slots. 
  \describe{
 	\item{access}{\code{.classname.slotname(object)}}
    \item{replace}{\code{.classname.slotname(object)<-}}
 	 }
 In this case classname is "registration" and slotname is for example "highres". So to access the "highres" slot of an r-object "myregistration"  (of class "registration") type \code{.registration.highres(myregistration)}.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\references{
See \url{http://www.fmrib.ox.ac.uk/fsl/} on how to obtain FSL and licences. See \url{http://www.fmrib.ox.ac.uk/fsl/data/atlas-descriptions.html} on information on the atlases.
}
\examples{
showClass("registration")
}
\keyword{classes}

