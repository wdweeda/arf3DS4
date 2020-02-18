\name{settings-class}
\Rdversion{1.1}
\docType{class}
\alias{settings-class}
\alias{.settings.avgdatFile,settings-method}
\alias{.settings.avgdatFile<-,settings-method}
\alias{.settings.avgDir,settings-method}
\alias{.settings.avgDir<-,settings-method}
\alias{.settings.avgtstatFile,settings-method}
\alias{.settings.avgtstatFile<-,settings-method}
\alias{.settings.avgWFile,settings-method}
\alias{.settings.avgWFile<-,settings-method}
\alias{.settings.betaDir,settings-method}
\alias{.settings.betaDir<-,settings-method}
\alias{.settings.conditionDir,settings-method}
\alias{.settings.conditionDir<-,settings-method}
\alias{.settings.conditionPrefix,settings-method}
\alias{.settings.conditionPrefix<-,settings-method}
\alias{.settings.dataDir,settings-method}
\alias{.settings.dataDir<-,settings-method}
\alias{.settings.dataRda,settings-method}
\alias{.settings.dataRda<-,settings-method}
\alias{.settings.derivativeFile,settings-method}
\alias{.settings.derivativeFile<-,settings-method}
\alias{.settings.expRda,settings-method}
\alias{.settings.expRda<-,settings-method}
\alias{.settings.funcDir,settings-method}
\alias{.settings.funcDir<-,settings-method}
\alias{.settings.funcRda,settings-method}
\alias{.settings.funcRda<-,settings-method}
\alias{.settings.logFile,settings-method}
\alias{.settings.logFile<-,settings-method}
\alias{.settings.lowresAvg,settings-method}
\alias{.settings.lowresAvg<-,settings-method}
\alias{.settings.lowresFile,settings-method}
\alias{.settings.lowresFile<-,settings-method}
\alias{.settings.modelDataFile,settings-method}
\alias{.settings.modelDataFile<-,settings-method}
\alias{.settings.modeldatDir,settings-method}
\alias{.settings.modeldatDir<-,settings-method}
\alias{.settings.modelDir,settings-method}
\alias{.settings.modelDir<-,settings-method}
\alias{.settings.modelnamesRda,settings-method}
\alias{.settings.modelnamesRda<-,settings-method}
\alias{.settings.modelPrefix,settings-method}
\alias{.settings.modelPrefix<-,settings-method}
\alias{.settings.modelRda,settings-method}
\alias{.settings.modelRda<-,settings-method}
\alias{.settings.optionsRda,settings-method}
\alias{.settings.optionsRda<-,settings-method}
\alias{.settings.regDir,settings-method}
\alias{.settings.regDir<-,settings-method}
\alias{.settings.regRda,settings-method}
\alias{.settings.regRda<-,settings-method}
\alias{.settings.residualFile,settings-method}
\alias{.settings.residualFile<-,settings-method}
\alias{.settings.startRda,settings-method}
\alias{.settings.startRda<-,settings-method}
\alias{.settings.statsDir,settings-method}
\alias{.settings.statsDir<-,settings-method}
\alias{.settings.statsRda,settings-method}
\alias{.settings.statsRda<-,settings-method}
\alias{.settings.subjectDir,settings-method}
\alias{.settings.subjectDir<-,settings-method}
\alias{.settings.subjectPrefix,settings-method}
\alias{.settings.subjectPrefix<-,settings-method}
\alias{.settings.version,settings-method}
\alias{.settings.version<-,settings-method}
\alias{.settings.weightFile,settings-method}
\alias{.settings.weightFile<-,settings-method}
\alias{.settings.weightsDir,settings-method}
\alias{.settings.weightsDir<-,settings-method}

\title{arf3DS4 "settings" class}
\description{The arf3DS4 class "settings" contains the names of all directories and files used by the arf3DS4 package. Default values are available for the "settings" class and adjusting the slot values is usually not necessary. An object of class "settings" is contained in the \code{\linkS4class{experiment}} class.}
\section{Objects}{
Objects can be created by calls of the form \code{new("settings", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{expRda}:}{Name of the .Rda file containing the experiment information.}
    \item{\code{optionsRda}:}{Name of the .Rda file containing the options for fitting an ARF model.}
    \item{\code{startRda}:}{Name of the .Rda file containing the startingvalues for fitting an ARF model.}
    \item{\code{dataRda}:}{Name of the .Rda file containing the fMRI data location and information.}
    \item{\code{modelRda}:}{Name of the .Rda file containing the ARF model information.}
    \item{\code{statsRda}:}{Name of the .Rda file containing the statistics information of a model (not yet implemented).}
    \item{\code{regRda}:}{Name of the .Rda file containing the registration information of fMRI data image.}
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
    \item{\code{funcDir}:}{Name of the functional directories (within subjects and conditions)/data.}
    \item{\code{betaDir}:}{Name of the beta directory (within conditions)/data.}
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
    \item{\code{lowresFile}:}{Name of the lowres structural image.}
    \item{\code{lowresAvg}:}{Name of the average lowres structural image.}
    \item{\code{logFile}:}{Name of the logFile (not yet implemented).}
    \item{\code{version}:}{Object of class "version" (see \code{\linkS4class{version}}) }
  }
}
\section{Methods}{
Standard arf3DS4 accessor and replacement functions can be used to access and replace slots.  
  \describe{
    \item{access}{\code{.classname.slotname(object)}}
    \item{replace}{\code{.classname.slotname(object)<-}}
	 }
In this case classname is "settings" and slotname is for example "avgdatfile". So to access the "avgdatfile" slot of an r-object "mysettings"  (of class "settings") type \code{.settings.avgdatfile(mysettings)}.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{experiment}}
}
\examples{
showClass("settings")
}
\keyword{classes}
