\name{experiment-class}
\Rdversion{1.1}
\docType{class}
\alias{experiment-class}
\alias{.experiment.avgdatFile,experiment-method}
\alias{.experiment.avgdatFile<-,experiment-method}
\alias{.experiment.avgDir,experiment-method}
\alias{.experiment.avgDir<-,experiment-method}
\alias{.experiment.avgtstatFile,experiment-method}
\alias{.experiment.avgtstatFile<-,experiment-method}
\alias{.experiment.avgWFile,experiment-method}
\alias{.experiment.avgWFile<-,experiment-method}
\alias{.experiment.betaDir,experiment-method}
\alias{.experiment.betaDir<-,experiment-method}
\alias{.experiment.condition.names,experiment-method}
\alias{.experiment.condition.names<-,experiment-method}
\alias{.experiment.condition.num,experiment-method}
\alias{.experiment.condition.num<-,experiment-method}
\alias{.experiment.conditionDir,experiment-method}
\alias{.experiment.conditionDir<-,experiment-method}
\alias{.experiment.conditionPrefix,experiment-method}
\alias{.experiment.conditionPrefix<-,experiment-method}
\alias{.experiment.dataDir,experiment-method}
\alias{.experiment.dataDir<-,experiment-method}
\alias{.experiment.dataRda,experiment-method}
\alias{.experiment.dataRda<-,experiment-method}
\alias{.experiment.derivativeFile,experiment-method}
\alias{.experiment.derivativeFile<-,experiment-method}
\alias{.experiment.expRda,experiment-method}
\alias{.experiment.expRda<-,experiment-method}
\alias{.experiment.funcDir,experiment-method}
\alias{.experiment.funcDir<-,experiment-method}
\alias{.experiment.funcRda,experiment-method}
\alias{.experiment.funcRda<-,experiment-method}
\alias{.experiment.logFile,experiment-method}
\alias{.experiment.logFile<-,experiment-method}
\alias{.experiment.lowresAvg,experiment-method}
\alias{.experiment.lowresAvg<-,experiment-method}
\alias{.experiment.lowresFile,experiment-method}
\alias{.experiment.lowresFile<-,experiment-method}
\alias{.experiment.modelDataFile,experiment-method}
\alias{.experiment.modelDataFile<-,experiment-method}
\alias{.experiment.modeldatDir,experiment-method}
\alias{.experiment.modeldatDir<-,experiment-method}
\alias{.experiment.modelDir,experiment-method}
\alias{.experiment.modelDir<-,experiment-method}
\alias{.experiment.modelnamesRda,experiment-method}
\alias{.experiment.modelnamesRda<-,experiment-method}
\alias{.experiment.modelPrefix,experiment-method}
\alias{.experiment.modelPrefix<-,experiment-method}
\alias{.experiment.modelRda,experiment-method}
\alias{.experiment.modelRda<-,experiment-method}
\alias{.experiment.name,experiment-method}
\alias{.experiment.name<-,experiment-method}
\alias{.experiment.optionsRda,experiment-method}
\alias{.experiment.optionsRda<-,experiment-method}
\alias{.experiment.path,experiment-method}
\alias{.experiment.path<-,experiment-method}
\alias{.experiment.regDir,experiment-method}
\alias{.experiment.regDir<-,experiment-method}
\alias{.experiment.regRda,experiment-method}
\alias{.experiment.regRda<-,experiment-method}
\alias{.experiment.residualFile,experiment-method}
\alias{.experiment.residualFile<-,experiment-method}
\alias{.experiment.startRda,experiment-method}
\alias{.experiment.startRda<-,experiment-method}
\alias{.experiment.statsDir,experiment-method}
\alias{.experiment.statsDir<-,experiment-method}
\alias{.experiment.statsRda,experiment-method}
\alias{.experiment.statsRda<-,experiment-method}
\alias{.experiment.subject.names,experiment-method}
\alias{.experiment.subject.names<-,experiment-method}
\alias{.experiment.subject.num,experiment-method}
\alias{.experiment.subject.num<-,experiment-method}
\alias{.experiment.subjectDir,experiment-method}
\alias{.experiment.subjectDir<-,experiment-method}
\alias{.experiment.subjectPrefix,experiment-method}
\alias{.experiment.subjectPrefix<-,experiment-method}
\alias{.experiment.version,experiment-method}
\alias{.experiment.version<-,experiment-method}
\alias{.experiment.weightFile,experiment-method}
\alias{.experiment.weightFile<-,experiment-method}
\alias{.experiment.weightsDir,experiment-method}
\alias{.experiment.weightsDir<-,experiment-method}
\alias{show,experiment-method}

\title{arf3DS4 "experiment" class}
\description{The "experiment" class is the backbone of the arf3DS4 package. The class contains information on the names of all relevant directories and files (these are inherited from an "settings" object (see \code{\linkS4class{settings}})). 
Before running an ARF analysis an object of class "experiment" must be loaded (using \code{\link{loadExp}}). Experiment directories can be created using \code{\link{makeExpDirs}}. To view information on an experiment call \code{\link{getExp}}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("experiment", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{path}:}{Path of the experiment. }
    \item{\code{name}:}{Name of the experiment. }
    \item{\code{subject.num}:}{Number of subjects. }
    \item{\code{subject.names}:}{Vector of \code{length(subject.num)} containing subject names. }
    \item{\code{condition.num}:}{Number of conditions.}
    \item{\code{condition.names}:}{Vector of \code{length(condition.num)} containing condition names. }
    \item{\code{expRda}:}{Name of the .Rda file containing the experiment information.}
    \item{\code{optionsRda}:}{Name of the .Rda file containing the options for fitting an ARF model.}
    \item{\code{startRda}:}{Name of the .Rda file containing the startingvalues for fitting an ARF model.}
    \item{\code{dataRda}:}{Name of the .Rda file containing the fMRI data location and information.}
    \item{\code{modelRda}:}{Name of the .Rda file containing the ARF model information.}
    \item{\code{statsRda}:}{Name of the .Rda file containing the statistics information of a model (not yet implemented).}
    \item{\code{regRda}:}{Name of the .Rda file containing the registration information of an fMRI data image.}
    \item{\code{funcRda}:}{Name of the .Rda file containing the functional information of an fMRI data image. }
    \item{\code{subjectPrefix}:}{Prefix for the subject names. Defaults to empty.}
    \item{\code{conditionPrefix}:}{Prefix for the condition names. Defaults to empty. }
    \item{\code{modelPrefix}:}{Prefix for the model names. Defaults to empty. }
    \item{\code{subjectDir}:}{Name of the subjects directory.}
    \item{\code{conditionDir}:}{Name of the conditions directory.}
    \item{\code{dataDir}:}{Name of the data directory.}
    \item{\code{weightsDir}:}{Name of the weights directory (within conditions/data).}
    \item{\code{avgDir}:}{Name of the average data directory (within conditions/data).}
    \item{\code{regDir}:}{Name of the registration directory (within conditions/data).}
    \item{\code{funcDir}:}{Name of the functional directories (within subjects and conditions/data).}
    \item{\code{betaDir}:}{Name of the beta directory (within conditions/data).}
    \item{\code{modelDir}:}{Name of the model directory.}
    \item{\code{statsDir}:}{Name of the stats directory (not yet implemented).}
    \item{\code{modeldatDir}:}{Name of the data directory (within models).}
    \item{\code{avgdatFile}:}{Name of the averaged data file (fMRI).}
    \item{\code{avgWFile}:}{Name of the averaged weight file (fMRI).}
    \item{\code{avgtstatFile}:}{Name of the averaged t-values file (fMRI).}
    \item{\code{modelDataFile}:}{Name of the model prediction file (fMRI).}
    \item{\code{modelnamesRda}:}{Name of the modelnames.Rda file containing the fitted models.}
    \item{\code{residualFile}:}{Name of the residuals file (binary).}
    \item{\code{derivativeFile}:}{Name of the derivatives file (binary).}
    \item{\code{weightFile}:}{Name of the weight file (binary).}
    \item{\code{lowresFile}:}{Name of the low resolution structural image.}
    \item{\code{lowresAvg}:}{Name of the average low resolution structural image.}
    \item{\code{logFile}:}{Name of the logFile (not yet implemented).}
    \item{\code{version}:}{Object of class "version" (see \code{\linkS4class{version}}) }
  }
}
\section{Extends}{
Class \code{"\linkS4class{settings}"}, directly.
}
\section{Methods}{
Standard arf3DS4 accessor and replacement functions can be used to access and replace slots.  
  \describe{
    \item{access}{\code{.classname.slotname(object)}}
    \item{replace}{\code{.classname.slotname(object)<-}}
	 }
In this case classname is "experiment" and slotname is for example "avgdatfile". So to access the "avgdatfile" slot of an r-object "myexperiment"  (of class "experiment") type \code{.experiment.avgdatfile(myexperiment)}.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{settings}}, \code{\link{makeExpDirs}}, \code{\link{loadExp}}, \code{\link{getExp}}
}

\examples{
showClass("experiment")
}
\keyword{classes}
