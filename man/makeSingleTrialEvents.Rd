\name{makeSingleTrialEvents}
\alias{makeSingleTrialEvents}

\title{
Estimate single-trial fMRI activity.
}
\description{
makeSingleTrialEvents calculates single-trial fMRI activity using a double gamma HRF.
}
\usage{
makeSingleTrialEvents(subject, condition, sefilename='single_events', 
 hrf.control=list(a1=6,a2=12,b1=0.9,b2=0.9,ce=0.35), experiment=NULL)
}

\arguments{
  \item{subject}{Name of the subject.}
  \item{condition}{Name of the condition.}
  \item{sefilename}{Filename of the volume where the single-trial data is stored (stores a 4D file with a volume for each event).}
  \item{hrf.control}{A list containing the parameters for the double gamma HRF.}
  \item{experiment}{An optional object of class \code{\linkS4class{experiment}}.}
}
\section{Warnings}{Before \code{makeSingleTrialEvents} can be executed, the data of each run within a condition must be linked to the raw time-series. This is
done via object of class \code{\linkS4class{functional}} by calling \code{\link{setFuncTimings}}.}

\details{
Files are saved in the \bold{data/functional} directory using the filename specified in \code{sefilename}. The function also returns the single-trial data. 
}
\value{
Returns a 4D nifti file (of class \code{\linkS4class{fmri.data}}) containing the single-trial data.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}, \code{\linkS4class{fmri.data}}, \code{\linkS4class{experiment}}, \code{\linkS4class{functional}}, \code{\link{setFuncTimings}}
}
\keyword{ models }
