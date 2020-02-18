\name{functional-class}
\Rdversion{1.1}
\docType{class}
\alias{functional-class}
\alias{.functional.filename,functional-method}
\alias{.functional.filename<-,functional-method}
\alias{.functional.fullpath,functional-method}
\alias{.functional.fullpath<-,functional-method}
\alias{.functional.functionaldata,functional-method}
\alias{.functional.functionaldata<-,functional-method}
\alias{.functional.linkedfile,functional-method}
\alias{.functional.linkedfile<-,functional-method}
\alias{.functional.linkedfiles,functional-method}
\alias{.functional.linkedfiles<-,functional-method}
\alias{.functional.timings,functional-method}
\alias{.functional.timings<-,functional-method}
\alias{.functional.version,functional-method}
\alias{.functional.version<-,functional-method}

\title{arf3DS4 "functional" class }
\description{The arf3DS4 "functional" class contains information on the functional volume where the fMRI beta-files were estimated on. It links the data of the different runs with the raw time-series files, primarily for connectivity analysis. 
To link runs to raw time-series call \code{\link{setFuncTimings}}.}
\section{Objects}{
Objects can be created by calls of the form \code{new("functional", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{fullpath}:}{Full path of the functional volume directory. }
    \item{\code{functionaldata}:}{Filename of the functional volume (raw time-series).}
    \item{\code{filename}:}{Filename of the functional.Rda file. }
    \item{\code{linkedfiles}:}{Path to linked beta-file.}
    \item{\code{timings}:}{Vector of stimulus timings.}
  	\item{\code{version}:}{Object of class "version" (see \code{\linkS4class{version}}) }
  }
}
\section{Methods}{
Standard arf3DS4 accessor and replacement functions can be used to access and replace slots.  
  \describe{
    \item{access}{\code{.classname.slotname(object)}}
    \item{replace}{\code{.classname.slotname(object)<-}}
	 }
In this case classname is "functional" and slotname is for example "timings". So to access the "timings" slot of an r-object "myfunctional"  (of class "functional") type \code{.functional.timings(myfunctional)}.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\examples{
showClass("functional")
}
\keyword{classes}
