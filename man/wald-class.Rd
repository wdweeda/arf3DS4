\name{wald-class}
\Rdversion{1.1}
\docType{class}
\alias{wald-class}
\alias{.wald.consts,wald-method}
\alias{.wald.consts<-,wald-method}
\alias{.wald.df1,wald-method}
\alias{.wald.df1<-,wald-method}
\alias{.wald.df2,wald-method}
\alias{.wald.df2<-,wald-method}
\alias{.wald.pvalues,wald-method}
\alias{.wald.pvalues<-,wald-method}
\alias{.wald.stats,wald-method}
\alias{.wald.stats<-,wald-method}

\title{arf3DS4 "wald" class}
\description{The arf3DS4 "wald" class contains information of the statistics on the spatial model parameters. It is used in the \code{\linkS4class{model}} class.}
\section{Objects}{
Objects can be created by calls of the form \code{new("wald", ...)}. 
}
\section{Slots}{
	 \describe{
    \item{\code{consts}:}{The #regions by #hypothesis matrix defining the constants used for testing.}
    \item{\code{stats}:}{The #regions by #hypothesis matrix of wald statistics.}
    \item{\code{df1}:}{Vector of the Degrees of Freedom of the model.}
    \item{\code{df2}:}{Vector of the Degrees of Freedom of the error. }
    \item{\code{pvalues}:}{The #regions by #hypothesis matrix of p-values.}
  }
}
\section{Methods}{
 Standard arf3DS4 accessor and replacement functions can be used to access and replace slots. 
  \describe{
 	\item{access}{\code{.classname.slotname(object)}}
    \item{replace}{\code{.classname.slotname(object)<-}}
	 }
 In this case classname is "wald" and slotname is for example "stats". So to access the "stats" slot of an r-object "mywald"  (of class "wald") type \code{.wald.stats(mywald)}.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}, \code{\link{wald}}
}
\examples{
showClass("wald")
}
\keyword{classes}

