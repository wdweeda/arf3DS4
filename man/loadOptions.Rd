\name{loadOptions}
\Rdversion{1.1}
\alias{loadOptions}
\alias{saveOptions}
\title{ Save/load ARF Options }
\description{
Load or save an "options" object belonging to an ARF model.
}
\usage{
loadOptions(arfmodel)
saveOptions(options,arfmodel)
}
\arguments{
  \item{arfmodel}{Object of class "model" (see \code{\linkS4class{model}})}
  \item{options}{Object of class "options" (see \code{\linkS4class{options}})}
  
}
\value{
loadOptions returns an object of class \code{\linkS4class{options}}.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\seealso{
	\code{\linkS4class{model}}, \code{\linkS4class{options}}
}
\keyword{ models }
