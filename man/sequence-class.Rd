\name{sequence-class}
\Rdversion{1.1}
\docType{class}
\alias{sequence-class}
\alias{.sequence.best,sequence-method}
\alias{.sequence.best<-,sequence-method}
\alias{.sequence.current,sequence-method}
\alias{.sequence.current<-,sequence-method}
\alias{.sequence.fit,sequence-method}
\alias{.sequence.fit<-,sequence-method}
\alias{.sequence.minimum,sequence-method}
\alias{.sequence.minimum<-,sequence-method}
\alias{.sequence.mnames,sequence-method}
\alias{.sequence.mnames<-,sequence-method}
\alias{.sequence.regions,sequence-method}
\alias{.sequence.regions<-,sequence-method}
\alias{.sequence.valid,sequence-method}
\alias{.sequence.valid<-,sequence-method}
\alias{show,sequence-method}

\title{arf3DS4 "sequence" class}
\description{The arf3DS4 "sequence" class contains information of a sequence of fitted models. It can be used to assess the relative fit between models and decide which model fits best.}
\section{Objects}{
Objects can be created by calls of the form \code{new("sequence", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{best}:}{Which model currently fits the best? }
    \item{\code{current}:}{The current model number in the sequence (not used).}
    \item{\code{regions}:}{Vector of number-of-regions to in the fitted sequence.}
    \item{\code{mnames}:}{Vector of modelnames of the sequence.}
    \item{\code{fit}:}{Vector of fit indices for each model.}
    \item{\code{minimum}:}{Vector of minima for each model.}
    \item{\code{valid}:}{Logical vector indicating valid models.}
  }
}
\section{Methods}{
 Standard arf3DS4 accessor and replacement functions can be used to access and replace slots. 
  \describe{
 	\item{access}{\code{.classname.slotname(object)}}
    \item{replace}{\code{.classname.slotname(object)<-}}
 	 }
 In this case classname is "sequence" and slotname is for example "fit". So to access the "fit" slot of an r-object "mysequence" (of class "sequence") type \code{.sequence.fit(mysequence)}.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\link{minBIC}}, \code{\link{fitRange}}
}
\examples{
showClass("sequence")
}
\keyword{classes}
