\name{data-class}
\Rdversion{1.1}
\docType{class}
\alias{data-class}
\alias{.data.avgdatfile,data-method}
\alias{.data.avgdatfile<-,data-method}
\alias{.data.avgtstatFile,data-method}
\alias{.data.avgtstatFile<-,data-method}
\alias{.data.avgWfile,data-method}
\alias{.data.avgWfile<-,data-method}
\alias{.data.betafiles,data-method}
\alias{.data.betafiles<-,data-method}
\alias{.data.fullpath,data-method}
\alias{.data.fullpath<-,data-method}
\alias{.data.funcDir,data-method}
\alias{.data.funcDir<-,data-method}
\alias{.data.funcRda,data-method}
\alias{.data.funcRda<-,data-method}
\alias{.data.mask,data-method}
\alias{.data.mask<-,data-method}
\alias{.data.n,data-method}
\alias{.data.n<-,data-method}
\alias{.data.name,data-method}
\alias{.data.name<-,data-method}
\alias{.data.regDir,data-method}
\alias{.data.regDir<-,data-method}
\alias{.data.regRda,data-method}
\alias{.data.regRda<-,data-method}
\alias{.data.ss,data-method}
\alias{.data.ss<-,data-method}
\alias{.data.runs,data-method}
\alias{.data.runs<-,data-method}
\alias{.data.dataHeader,data-method}
\alias{.data.dataHeader<-,data-method}
\alias{.data.version,data-method}
\alias{.data.version<-,data-method}
\alias{.data.weightfiles,data-method}
\alias{.data.weightfiles<-,data-method}
\alias{show,data-method}

\title{arf3DS4 "data" class}
\description{The "data" class contains all information on the fMRI data of one condition of a subject. Specifically it contains filenames of the datafiles (with full paths) 
and some basic information of the data, the number of runs, and additional files used by ARF. Each instance of a \code{\linkS4class{model}} class object inherits a "data" object. 
The "data" class also uses some slots of the \code{\linkS4class{experiment}} class.}
\section{Objects}{
Objects can be created by calls of the form \code{new("data", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{name}:}{Name of the data object, usually 'subjectname - conditionname'. }
    \item{\code{fullpath}:}{Full path to the data directory.}
    \item{\code{betafiles}:}{Vector of filenames (full paths) of the beta files (containing beta estimates).}
    \item{\code{weightfiles}:}{Vector of filenames (full path) of the weight files (containing variances of the beta estimates).}
    \item{\code{avgdatfile}:}{Filename of the averaged beta file.}
    \item{\code{avgWfile}:}{Filename of the averaged weight file.}
    \item{\code{avgtstatFile}:}{Filename of the average t-statistics (=avgbetafile/sqrt(avgweightfile)) file. }
    \item{\code{n}:}{Number of 'brain' voxels (excluding voxels not masked).}
    \item{\code{mask}:}{Vector containing the mask for the data.}
    \item{\code{ss}:}{Sums-of-squares of the data.}
    \item{\code{regDir}:}{Full path to the registration directory.}
    \item{\code{regRda}:}{Name of the registration.Rda file.}
    \item{\code{funcDir}:}{Full path to the functional directory.}
    \item{\code{funcRda}:}{Name of the functional.Rda file.}
    \item{\code{dataHeader}:}{Nifti header information from the avgdataFile.}
    \item{\code{runs}:}{Number of runs (= number of beta files).}
    \item{\code{version}:}{Object of class "version" (see \code{\linkS4class{version}}).}
  }
}
\details{
The "data" class contains information on the file-locations of the fMRI data. The fMRI data (outcomes of a GLM analysis) must be available in the \bold{/beta} and \bold{/weights} directories of each condition for each subject. 
The data can be in the form of beta-values with associated squared standard errors or t-values with the standard errors set to 1. In the latter case the only files that have to be copied are the t-values to the \bold{/beta} directory,
the weight-files are then created automatically. The average files are created after a call to \code{\link{createAverages}}.
}
\section{Methods}{
Standard arf3DS4 accessor and replacement functions can be used to access and replace slots.  
  \describe{
    \item{access}{\code{.classname.slotname(object)}}
    \item{replace}{\code{.classname.slotname(object)<-}}
	 }
In this case classname is "data" and slotname is for example "avgstatFile". So to access the "avgtstatFile" slot of an r-object "mydata"  (of class "data") type \code{.settings.avgtstatFile(mydata)}.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}, \code{\linkS4class{experiment}}
}
\examples{
showClass("data")
}
\keyword{classes}
