\name{plot.vscc}
\alias{plot.vscc}
\title{Plotting for VSCC Objects}
\description{Dedicated plot function for objects of class vscc.}
\usage{
\method{plot}{vscc}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class vscc.
}
\item{\dots}{
Further arguments to be passed on
}
}
\details{
Provides a scatterplot matrix of the selected variables with colours corresponding to each group.
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
X<-banknote[,-1]
bankrun <- vscc(X)
plot(bankrun)
}
%\keyword{ ~kwd1 }
