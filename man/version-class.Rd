\name{version-class}
\Rdversion{1.1}
\docType{class}
\alias{version-class}
\alias{.version.build,version-method}
\alias{.version.build<-,version-method}
\alias{.version.svnrev,version-method}
\alias{.version.svnrev<-,version-method}
\alias{.version.update,version-method}
\alias{.version.update<-,version-method}
\alias{.version.version,version-method}
\alias{.version.version<-,version-method}
\alias{show,version-method}

\title{arf3DS4 "version" class}
\description{The "version" class contains information on the version of the arf3DS4 objects. Each arf3DS4 object contains this class. }
\section{Objects}{
Objects can be created by calls of the form \code{new("version", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{version}:}{Version number of the arf3DS4 package.}
    \item{\code{build}:}{Build number of the arf3DS4 package.}
    \item{\code{update}:}{Update number of the arf3DS4 package.}
    \item{\code{svnrev}:}{SVN revision on r-forge.}
  }
}
\section{Methods}{
Standard arf3DS4 accessor and replacement functions can be used to access and replace slots.  
  \describe{
    \item{access}{\code{.classname.slotname(object)}}
    \item{replace}{\code{.classname.slotname(object)<-}}
	 }
In this case classname is "version" and slotname is for example "snvrev". So to access the "svnrev" slot of an r-object "currentversion"  (of class "version") type \code{.version.svnrev(currentversion)}.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\note{ The "version" class is mainly used for back-compatibility issues.}
\examples{
showClass("version")
}
\keyword{classes}
