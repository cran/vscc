\name{summary.vscc}
\alias{summary.vscc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary for VSCC Objects
}
\description{
Dedicated summary function for objects of class vscc
}
\usage{
\method{summary}{vscc}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class vscc
}
\item{\dots}{
Additional arguments to be passed
}
}
\value{
No return value.
}
\author{
Jeffrey L. Andrews
}


\seealso{
\code{\link{vscc}}
}
\examples{
require("mclust")
data(banknote)
summary(vscc(banknote[,-1]))
}
%\keyword{ ~kwd1 }

