\name{mnames-class}
\Rdversion{1.1}
\docType{class}
\alias{mnames-class}
\alias{.mnames.condition,mnames-method}
\alias{.mnames.condition<-,mnames-method}
\alias{.mnames.experiment,mnames-method}
\alias{.mnames.experiment<-,mnames-method}
\alias{.mnames.mnames,mnames-method}
\alias{.mnames.mnames<-,mnames-method}
\alias{.mnames.subject,mnames-method}
\alias{.mnames.subject<-,mnames-method}
\alias{show,mnames-method}

\title{arf3DS4 "mnames" class}
\description{The arf3DS4 "mnames" class contains information on the models in a subject/condition. It can be used to quickly load models without having to specify the experiment, subject and condition each time.}
\section{Objects}{
Objects can be created by calls of the form \code{new("mnames", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{experiment}:}{Object of class "experiment" (see \code{\linkS4class{experiment}}).}
    \item{\code{subject}:}{Subject name.}
    \item{\code{condition}:}{Condition name.}
    \item{\code{mnames}:}{Vector of modelnames within a subject/condition.}
  }
}
\section{Methods}{
 Standard arf3DS4 accessor and replacement functions can be used to access and replace slots. 
  \describe{
 	\item{access}{\code{.classname.slotname(object)}}
    \item{replace}{\code{.classname.slotname(object)<-}}
    \item{show}{\code{signature(object = "mnames")}: ... }
	 }
 In this case classname is "mnames" and slotname is for example "subject". So to access the "subject" slot of an r-object "mymnames"  (of class "mnames") type \code{.mnames.subject(mymnames)}.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{experiment}}
}
\examples{
showClass("mnames")
}
\keyword{classes}
