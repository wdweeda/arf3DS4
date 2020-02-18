\name{options-class}
\Rdversion{1.1}
\docType{class}
\alias{options-class}
\alias{.options.adjust.n,options-method}
\alias{.options.adjust.n<-,options-method}
\alias{.options.chk.method,options-method}
\alias{.options.chk.method<-,options-method}
\alias{.options.chk.range,options-method}
\alias{.options.chk.range<-,options-method}
\alias{.options.min.analyticalgrad,options-method}
\alias{.options.min.analyticalgrad<-,options-method}
\alias{.options.min.iterlim,options-method}
\alias{.options.min.iterlim<-,options-method}
\alias{.options.min.boundlim,options-method}
\alias{.options.min.boundlim<-,options-method}
\alias{.options.min.routine,options-method}
\alias{.options.min.routine<-,options-method}
\alias{.options.nlm.gradtol,options-method}
\alias{.options.nlm.gradtol<-,options-method}
\alias{.options.nlm.steptol,options-method}
\alias{.options.nlm.steptol<-,options-method}
\alias{.options.opt.lower,options-method}
\alias{.options.opt.lower<-,options-method}
\alias{.options.opt.method,options-method}
\alias{.options.opt.method<-,options-method}
\alias{.options.opt.upper,options-method}
\alias{.options.opt.upper<-,options-method}
\alias{.options.output.mode,options-method}
\alias{.options.output.mode<-,options-method}
\alias{.options.start.maxfac,options-method}
\alias{.options.start.maxfac<-,options-method}
\alias{.options.start.method,options-method}
\alias{.options.start.method<-,options-method}
\alias{.options.start.vector,options-method}
\alias{.options.start.vector<-,options-method}
\alias{.options.sw.type,options-method}
\alias{.options.sw.type<-,options-method}
\alias{.options.version,options-method}
\alias{.options.version<-,options-method}
\alias{show,options-method}


\title{arf3DS4 "options" class  }
\description{The arf3DS4 "options" class contains the options of the minimization procedure (including starting value calculation, and (co)variance estimation). Mostly used by ARF "model" objects (see \code{\linkS4class{model}}).}
\section{Objects}{
Objects can be created by calls of the form \code{new("options", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{nlm.gradtol}:}{<nlm> Gradient tolerance (not used). }
    \item{\code{nlm.steptol}:}{<nlm> Step tolerance (not used). }
    \item{\code{opt.method}:}{<optim> Optim method (L-BFGS-B). }
    \item{\code{opt.lower}:}{<optim> Lowerbound of parameters (for L-BFGS-B).}
    \item{\code{opt.upper}:}{<optim> Upperbound of parameters (for L-BFGS-B).}
    \item{\code{min.analyticalgrad}:}{Use analytical gradient in minimization? }
    \item{\code{min.iterlim}:}{Iteration limit.}
    \item{\code{min.boundlim}:}{Persistent Boundary iteration limit.}
    \item{\code{min.routine}:}{Minimization routine (optim). }
    \item{\code{start.method}:}{Method of startingvalue calculation (\code{'use'},\code{'load'},or \code{'rect'})}
    \item{\code{start.maxfac}:}{Width of smoothing kernel used in starting value calculation.}
    \item{\code{start.vector}:}{Starting value template (used for default values).}
    \item{\code{chk.method}:}{Check method for parameters of minimized solution.}
    \item{\code{chk.range}:}{Range of values to check parameters of minimized solution.}
    \item{\code{sw.type}:}{Use diagonal residual matrix or full residual matrix in sandwich (co)variance calculation. }
    \item{\code{output.mode}:}{Vector with output types (\code{'none'} or \code{'progress'}).}
    \item{\code{version}:}{Object of class "version" (see \code{\linkS4class{version}})}
  }
}
\section{Methods}{
 Standard arf3DS4 accessor and replacement functions can be used to access and replace slots. 
  \describe{
 	\item{access}{\code{.classname.slotname(object)}}
    \item{replace}{\code{.classname.slotname(object)<-}}
 	 }
 In this case classname is "options" and slotname is for example "sw.type". So to access the "sw.type" slot of an r-object "myoptions"  (of class "options") type \code{.options.sw.type(myoptions)}.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}
}
\examples{
showClass("options")
}
\keyword{classes}

