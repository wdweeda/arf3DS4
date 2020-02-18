\name{newModel}
\Rdversion{1.1}
\alias{newModel}
\title{Create a new ARF model}
\description{
newModel creates a new model of with a given number of regions in the spatial model within a subject/condition, using the options specified. It also creates a new directory in the \bold{model} directory of the condition.
}
\usage{
newModel(modelname = "defaultmodel", regions = 1, subject = "", 
 condition = "", type = c("gauss", "simple"), options = new("options"), 
 overwrite = T, experiment = NULL)
}
\arguments{
  \item{modelname}{Name of the model.}
  \item{regions}{Number of regions in the spatial model.}
  \item{subject}{Name of the subject.}
  \item{condition}{Name of the condition.}
  \item{type}{Type of model to fit (full Gaussian (\code{'gauss'}), simple Gaussian (\code{'simple'}).}
  \item{options}{Options object (\code{\linkS4class{options}}).}
  \item{overwrite}{Overwrite the model if a model with the same name alreay exists?}
  \item{experiment}{Experiment object (as default uses \code{.experiment}).}
}
\details{
Calls to newModel save an instance of a file named model.Rda in the appropriate directory and return the model object.
}
\value{
Returns an object of class "model" with the appropriate slots.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
 \code{\linkS4class{model}}, \code{\linkS4class{options}}, \code{\linkS4class{experiment}}
 }
\keyword{ models }

