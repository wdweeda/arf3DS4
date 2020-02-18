\name{arf3DS4-package}
\Rdversion{1.1}
\alias{arf3DS4-package}
\alias{arf3DS4}
\docType{package}
\title{
Activated Region Fitting fMRI data analysis
}
\description{
Activated Region Fitting (ARF) is an analysis method for fMRI data. The ARF method uses Gaussian shaped functions to model active brain regions. 
An entire fMRI volume can be described by multiple parameterized Gaussian regions of activation. This parameterization allows for hypotheses on the
location of an active region, the spatial extent of an active region, and the amplitude of an active region. ARF can be invoked after standard GLM analysis, 
and uses the beta-values or t-values from this analysis. In addition ARF can be used to estimate functional connectivity, for this the raw time-series
(on which the GLM was performed) are also needed. 
The arf3DS4 package can read and write fMRI data in Nifti format (\code{\linkS4class{fmri.data}}). It is build to be compatible with fMRI analysis packages (e.g. FSL or SPM). The package 
uses a predefined directory- and file-structure to store its data. All objects (S4-class) and fMRI data files are stored in this directory- and file-structure.
}
\section{Quick start}{
The main directory- and file-structure (termed an experiment) is handled by an object of class \code{\linkS4class{experiment}}. To make an empty experiment structure
call \code{\link{makeExpDirs}}. This will create a directory- and file-structure to hold the fMRI data (which will have to be copied manually to the appropriate directories, see \code{\linkS4class{data}} for an overview of the data-directories ).
Once this is done an experiment is loaded by calling \code{\link{loadExp}}. 
To fit the ARF models (class \code{\linkS4class{model}}) to the data, first create a model using \code{\link{newModel}}, adjust the options of the model via the \code{\linkS4class{options}}-object using \code{\link{loadOptions}} and \code{\link{saveOptions}},
and finally fit the model using \code{\link{fitModel}}. To perform hypothesis tests call \code{\link{varcov}} and \code{\link{wald}} to calculate the parameter (co)variance matrix and Wald statistics respectively.
In addition to the hypothesis test connectivity between brain regions can be estimated (call \code{\link{fitConnectivity}} to perform this analysis). This function requires that the 
raw time-series of the fMRI data are available and that single-trial estimates (call \code{\link{makeSingleTrialEvents}}) are calculated.
}
\section{Example data}{
The package comes with an example dataset (load it using \code{data("arf-example-data")}). By then calling \code{\link{makeExample}}, a directory-structure containing
the example data is then created and loaded. Call \code{\link{getExp}} to view the data within the experiment.
}
\references{
\tabular{l}{
Wouter D. Weeda, Lourens J. Waldorp, Ingrid Christoffels and Hilde M. Huizenga (2009). Activated Region Fitting: A Robust High-Power Method for fMRI Analysis Using Parameterized Regions of Activation. Human Brain Mapping, 30(8), 2595-2605.\cr
Wouter D. Weeda, Lourens J. Waldorp, Raoul P.P.P. Grasman, Simon van Gaal and Hilde M. Huizenga (2011). Functional Connectivity Analysis of fMRI Data Using Parameterized Regions of Interest. NeuroImage, 54(1), 410-416.\cr
Wouter D. Weeda, Frank de Vos, Lourens J. Waldorp, Raoul P.P.P. Grasman, Hilde M. Huizenga (2011). arf3DS4: An Integrated Framework for Localization and Connectivity Analysis of fMRI Data. Journal of Statistical Software, 44(14), 1-33. URL \url{http://www.jstatsoft.org/v44/i14/} \cr
}
}

\details{
\tabular{ll}{
Package: \tab arf3DS4\cr
Type: \tab Package\cr
Version: \tab 2.5-4\cr
Date: \tab 2011-10-28\cr
Depends: \tab R (>= 2.12.0), methods, graphics, tcltk\cr
Repository/R-Forge/Revision: \tab 203\cr
Publication: \tab 2011-10-28 12:38:50\cr
Packaged: \tab 2011-10-28 12:38:50; rforge\cr
Built: \tab R 2.12.0; universal-apple-darwin9.8.0; 2011-10-28 12:38:50; unix\cr
LazyLoad: \tab yes\cr
}
}
\author{
Wouter D. Weeda <w.d.weeda@gmail.com>
}
\keyword{ package }
