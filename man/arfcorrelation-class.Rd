\name{arfcorrelation-class}
\Rdversion{1.1}
\docType{class}
\alias{arfcorrelation-class}
\alias{.arfcorrelation.corr,arfcorrelation-method}
\alias{.arfcorrelation.corr.pval,arfcorrelation-method}
\alias{.arfcorrelation.corr.pval<-,arfcorrelation-method}
\alias{.arfcorrelation.corr<-,arfcorrelation-method}
\alias{.arfcorrelation.num.corr,arfcorrelation-method}
\alias{.arfcorrelation.num.corr<-,arfcorrelation-method}
\alias{.arfcorrelation.pacorr,arfcorrelation-method}
\alias{.arfcorrelation.pacorr.pval,arfcorrelation-method}
\alias{.arfcorrelation.pacorr.pval<-,arfcorrelation-method}
\alias{.arfcorrelation.pacorr<-,arfcorrelation-method}
\alias{.arfcorrelation.timebyreg,arfcorrelation-method}
\alias{.arfcorrelation.timebyreg<-,arfcorrelation-method}

\title{Class "arfcorrelation"}
\description{ Object of class "arfcorrelation" hold information on correlations between regions of the spatial model. The functions \code{\link{fitConnectivity}} and \code{\link{roiConnectivity}} have as output an object of class "arfcorrelation".}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("arfcorrelation", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{timebyreg}:}{Trial-by-trial amplitude estimates per region (time x nreg)}
    \item{\code{corr}:}{Correlation matrix.}
    \item{\code{corr.pval}:}{P-values of the correlations}
    \item{\code{pacorr}:}{Partial correlation matrix (corrected for all other regions).}
    \item{\code{pacorr.pval}:}{P-values of partial correlations}
    \item{\code{num.corr}:}{Number of correlations tested.}
  }
}

\section{Methods}{
Standard arf3DS4 accessor and replacement functions can be used to access and replace slots.  
  \describe{
    \item{access}{\code{.classname.slotname(object)}}
    \item{replace}{\code{.classname.slotname(object)<-}}
	 }
In this case classname is "arfcorrelation" and slotname is for example "correlation". So to access the "correlation" slot of an r-object "arfcor"  (of class "arfcorrelation") type \code{.arfcorrelation.correlation(arfcor)}.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\examples{
showClass("arfcorrelation")
}
\keyword{classes}

