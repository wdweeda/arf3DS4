\name{differenceCor}
\alias{differenceCor}

\title{Calculate if correlations differ.}
\description{Calculate which correlations of two correlation matrices are different.}
\usage{
differenceCor(c1, c2, n1, n2 = n1)
}
\arguments{
  \item{c1}{Correlation matrix of sample 1.}
  \item{c2}{Correlation matrix of sample 2.}
  \item{n1}{Size of sample 1.}
  \item{n2}{Size of sample 2.}
}
\value{
  \item{z1}{Matrix of Fischer Transformed r's.}
  \item{z2}{Matrix of Fischer Transformed r's.}
  \item{dif}{Difference scores (z1-z2).}
  \item{z}{Z-value of the difference.}
  \item{pval }{P-value of z-diff.}
}
\author{ Wouter D. Weeda - w.d.weeda@gmail.com }
\keyword{ models }
