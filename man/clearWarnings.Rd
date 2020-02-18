\name{clearWarnings}
\Rdversion{1.1}
\alias{clearWarnings}
\title{ Clear warnings of a model object.}
\description{
Clearwarnings clears all warnings of a model object, sets the valid slot to \code{TRUE} and saves the model object. It asks for confirmation.
}
\usage{
clearWarnings(arfmodel, resetValid = T)
}
\arguments{
  \item{arfmodel}{Object of class "model" (see \code{\linkS4class{model}})}
  \item{resetValid}{If \code{TRUE} set the \code{@valid} slot to \code{TRUE}, else the \code{@valid} slot is left unchanged.}
}
\details{
Warnings are only cleared when starting the minimization procedure on a model object. After minimization all warnings are added to the \code{warnings} slot.
}
\value{
None
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}
}
\keyword{ models }

