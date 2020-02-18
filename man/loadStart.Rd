\name{loadStart}
\Rdversion{1.1}
\alias{loadStart}
\alias{saveStart}
\title{ Save/load ARF starting values }
\description{
Load the starting values of an ARF model from the .Rda file
}
\usage{
loadStart(arfmodel)
saveStart(startval,arfmodel)
}
\arguments{
  \item{arfmodel}{Object of class "model" (see \code{\linkS4class{model}})}
  \item{startval}{Vector of starting values.}
}
\value{
loadStart returns a vector with starting values.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}
}
\keyword{ models }