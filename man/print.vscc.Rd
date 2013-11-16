\name{print.vscc}
\alias{print.vscc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Printing for vscc
}
\description{
Dedicated print function for objects of class vscc.
}
\usage{
\method{print}{vscc}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class vscc
}
  \item{\dots}{
Further arguments to be passed on
}
}
\details{
Same as summary.
}
\author{
Jeffrey L. Andrews
}

\seealso{
\code{\link{summary.vscc}}, \code{\link{vscc}}
}
\examples{
require("mclust")
data(banknote)
vscc(banknote[,-1])
}
\keyword{ ~kwd1 }
