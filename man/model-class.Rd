\name{model-class}
\Rdversion{1.1}
\docType{class}
\alias{model-class}
\alias{.model.avgdatfile,model-method}
\alias{.model.avgdatfile<-,model-method}
\alias{.model.avgtstatFile,model-method}
\alias{.model.avgtstatFile<-,model-method}
\alias{.model.avgWfile,model-method}
\alias{.model.avgWfile<-,model-method}
\alias{.model.betafiles,model-method}
\alias{.model.betafiles<-,model-method}
\alias{.model.convergence,model-method}
\alias{.model.convergence<-,model-method}
\alias{.model.derivativeFile,model-method}
\alias{.model.derivativeFile<-,model-method}
\alias{.model.estimates,model-method}
\alias{.model.estimates<-,model-method}
\alias{.model.fit,model-method}
\alias{.model.fit<-,model-method}
\alias{.model.fullmodelDataFile,model-method}
\alias{.model.fullmodelDataFile<-,model-method}
\alias{.model.fullpath,model-method}
\alias{.model.fullpath<-,model-method}
\alias{.model.funcDir,model-method}
\alias{.model.funcDir<-,model-method}
\alias{.model.funcRda,model-method}
\alias{.model.funcRda<-,model-method}
\alias{.model.gradient,model-method}
\alias{.model.gradient<-,model-method}
\alias{.model.hessian,model-method}
\alias{.model.hessian<-,model-method}
\alias{.model.iterates,model-method}
\alias{.model.iterates<-,model-method}
\alias{.model.logFile,model-method}
\alias{.model.logFile<-,model-method}
\alias{.model.mask,model-method}
\alias{.model.mask<-,model-method}
\alias{.model.minimum,model-method}
\alias{.model.minimum<-,model-method}
\alias{.model.modelDataFile,model-method}
\alias{.model.modelDataFile<-,model-method}
\alias{.model.modeldatapath,model-method}
\alias{.model.modeldatapath<-,model-method}
\alias{.model.modelFile,model-method}
\alias{.model.modelFile<-,model-method}
\alias{.model.modelname,model-method}
\alias{.model.modelname<-,model-method}
\alias{.model.modelpath,model-method}
\alias{.model.modelpath<-,model-method}
\alias{.model.modeltype,model-method}
\alias{.model.modeltype<-,model-method}
\alias{.model.n,model-method}
\alias{.model.n<-,model-method}
\alias{.model.name,model-method}
\alias{.model.name<-,model-method}
\alias{.model.optionsFile,model-method}
\alias{.model.optionsFile<-,model-method}
\alias{.model.params,model-method}
\alias{.model.params<-,model-method}
\alias{.model.proctime,model-method}
\alias{.model.proctime<-,model-method}
\alias{.model.regDir,model-method}
\alias{.model.regDir<-,model-method}
\alias{.model.regions,model-method}
\alias{.model.regions<-,model-method}
\alias{.model.regRda,model-method}
\alias{.model.regRda<-,model-method}
\alias{.model.residualFile,model-method}
\alias{.model.residualFile<-,model-method}
\alias{.model.sandwichmethod,model-method}
\alias{.model.sandwichmethod<-,model-method}
\alias{.model.ss,model-method}
\alias{.model.ss<-,model-method}
\alias{.model.startFile,model-method}
\alias{.model.startFile<-,model-method}
\alias{.model.startval,model-method}
\alias{.model.startval<-,model-method}
\alias{.model.runs,model-method}
\alias{.model.runs<-,model-method}
\alias{.model.valid,model-method}
\alias{.model.valid<-,model-method}
\alias{.model.varcov,model-method}
\alias{.model.varcov<-,model-method}
\alias{.model.version,model-method}
\alias{.model.version<-,model-method}
\alias{.model.wald,model-method}
\alias{.model.wald<-,model-method}
\alias{.model.warnings,model-method}
\alias{.model.warnings<-,model-method}
\alias{.model.weightFile,model-method}
\alias{.model.weightFile<-,model-method}
\alias{.model.weightfiles,model-method}
\alias{.model.weightfiles<-,model-method}
\alias{.model.dataHeader,model-method}
\alias{.model.dataHeader<-,model-method}
\alias{plot,model,missing-method}
\alias{show,model-method}

\title{arf3DS4 "model" class}
\description{The arf3DS4 model class contains all information on a (fitted) model. It extends the "data" class directly, since it uses the files specified there to fit the model. 
The "model" class contains information on file locations and filenames. More importantly it contains information on the fitted model: parameter estimates, 
parameter (co)variance matrix, statistical test information, fit indices, convergence information, warnings, etc.}
\section{Objects}{
Objects can be created by calls of the form \code{new("model", ...)}.
 
}
\section{Slots}{
	 \describe{
    \item{\code{modelname}:}{Name of the (fitted) model.}
    \item{\code{modelpath}:}{Full path to the model directory.}
    \item{\code{modeldatapath}:}{Full path to the model/data directory. }
    \item{\code{residualFile}:}{Name of the file containing model residuals. }
    \item{\code{derivativeFile}:}{Name of the file containing first order derivatives. }
    \item{\code{weightFile}:}{Name of the binary file containing weights.}
    \item{\code{modelDataFile}:}{Name of the file containing model predictions. }
    \item{\code{fullmodelDataFile}:}{Full name of the file containing model predictions.}
    \item{\code{modelFile}:}{Name of the model.Rda file containing the model object.}
    \item{\code{optionsFile}:}{Name of the options.Rda file (used to specify options, see \code{\linkS4class{options}}).}
    \item{\code{startFile}:}{Name of the start.Rda file.}
    \item{\code{logFile}:}{Name of the log.Rda file. }
    \item{\code{convergence}:}{Convergence information.}
    \item{\code{iterates}:}{Number of iterates. }
    \item{\code{minimum}:}{Value of the objective function at the minimum. }
    \item{\code{estimates}:}{Parameter Estimates. }
    \item{\code{gradient}:}{Gradient of the parameters at the minimum. }
    \item{\code{hessian}:}{Hessian matrix of parameters.}
    \item{\code{params}:}{Number of parameters. }
    \item{\code{modeltype}:}{Type of model fitted.}
    \item{\code{sandwichmethod}:}{Method used to calculate Sandwich estimates.}
    \item{\code{varcov}:}{Variance/covariance matrix of the parameters. }
    \item{\code{warnings}:}{Warnings encountered during processing. }
    \item{\code{fit}:}{Fit of the model (BIC and RMSEA). }
    \item{\code{wald}:}{Wald class object (containing model statistics, see \code{\linkS4class{wald}}).}
    \item{\code{regions}:}{Number of regions in the spatial model.}
    \item{\code{startval}:}{Starting values of the minimization procedure. }
    \item{\code{proctime}:}{Processing time of minimization and sandwich procedures.}
    \item{\code{valid}:}{Is the model valid?}
    \item{\code{name}:}{Name of the data object, usually 'subjectname - conditionname'. }
    \item{\code{fullpath}:}{Full path to the data directory.}
    \item{\code{betafiles}:}{Vector of filenames (full paths) of the beta files (containing beta estimates).}
    \item{\code{weightfiles}:}{Vector of filenames (full paths) of the weight files (containing variances of the beta estimates).}
    \item{\code{avgdatfile}:}{Filename of the averaged beta file.}
    \item{\code{avgWfile}:}{Filename of the averaged weight file.}
    \item{\code{avgtstatFile}:}{Filename of the average t-statistics (=avgbetafile/sqrt(avgweightfile)) file. }
    \item{\code{n}:}{Number of 'brain' voxels (i.e. the number of voxels excluding those where \code{@mask} is zero).}
    \item{\code{mask}:}{Vector containing the mask for the data.}
    \item{\code{ss}:}{Sums-of-squares of the data.}
    \item{\code{regDir}:}{Fullpath to the registration directory.}
    \item{\code{regRda}:}{Name of the registration.Rda file.}
    \item{\code{funcDir}:}{Fullpath to the functional directory.}
    \item{\code{funcRda}:}{Name of the functional.Rda file.}
    \item{\code{dataHeader}:}{Nifti header information from the avgdataFile.}
    \item{\code{runs}:}{Number of runs (= number of beta files).}
    \item{\code{version}:}{Object of class "version" (see \code{\linkS4class{version}})}
  }
}
\section{Extends}{
Class \code{\linkS4class{data}}, directly.
}
\section{Methods}{
 Standard arf3DS4 accessor and replacement functions can be used to access and replace slots. 
  \describe{
 	\item{access}{\code{.classname.slotname(object)}}
    \item{replace}{\code{.classname.slotname(object)<-}}
    \item{plot}{\code{signature(x = "model", y = "missing")}: ... }
    \item{show}{\code{signature(object = "model")}: ... }
	 }
 In this case classname is "model" and slotname is for example "estimates". So to access the "estimates" slot of an r-object "mymodel"  (of class "model") type \code{.model.estimates(mymodel)}.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{data}}, \code{\linkS4class{options}}, \code{\linkS4class{wald}}
}
\examples{
showClass("model")
}
\keyword{classes}
