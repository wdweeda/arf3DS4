\name{loadRda}
\Rdversion{1.1}
\alias{loadRda}
\title{Load an .Rda file into a specified object.}
\description{
loadRda loads a standard .Rda file, and returns this object (instead of loading it directly in the R workspace).
}
\usage{
loadRda(file)
}
\arguments{
  \item{file}{Filename}
}
\details{
loadRda assumes that only one object was saved. It cannot load an entire workspace.
}
\value{
Returns the object that was saved in the .Rda file.
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\keyword{ file }

