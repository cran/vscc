\name{print.vscc}
\alias{print.vscc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Printing for VSCC
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
\value{
No return value.
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
X<-banknote[,-1]
vscc(X)
}
%\keyword{ ~kwd1 }
