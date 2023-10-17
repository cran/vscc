\name{vscc}
\alias{vscc}
\title{Variable Selection for Clustering and Classification}
\description{
Performs variable selection under a clustering or classification framework. Automated implementation using model-based clustering is based on \code{teigen} version 2.0 and \code{mclust} version 4.0; issues *may* arise when using different versions.
}
\usage{
vscc(x, G=1:9, automate = "mclust", initial = NULL, initunc=NULL, train = NULL,
    forcereduction = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data frame or matrix to perform variable selection on
}
 \item{G}{
Vector for the number of groups to consider during initialization and/or post-selection analysis. Default is 1-9. 
}
  \item{automate}{
Character string (\code{"teigen"}, \code{"mclust"} (default), or NULL only) indicating which mixture model family to implement as initialization and/or post-selection analysis. If NULL, the function assumes manual operation of the algorithm (meaning an initial clustering vector must be given, and no post-selection analysis is performed).
}
  \item{initial}{
Optional vector giving the initial clustering. 
}
  \item{initunc}{
Optional scalar indicating the total uncertainty of the initial clustering solution. Only used when \code{initial} is non-null. 
}
  \item{train}{
Optional vector of training data (for classification framework).
}
  \item{forcereduction}{
Logical indicating if the full data set should be considered (FALSE) when selecting the `best' variable subset via total model uncertainty. Not used if \code{automate=NULL}.
}
}
\value{
\item{selected}{
A list containing the subsets of variables selected for each relation. Each set is numbered according to the number in the exponential of the relationship. For instance, \code{vscc_object$selected[[3]]} corresponds to the variable subset selected by the cubic relationship.
}
\item{family}{
The family used as initialization and/or post selection. (Same as user input \code{automate}, and can be \code{NULL}).
}
\item{wss}{
The within-group variance associated with each variable from the full data set.
}

The remaining values are provided as long as \code{automate} is not \code{NULL}:

\item{topselected}{
The best variable subset according to the total model uncertainty.
}
\item{initialrun}{
Results from the initialization; an object of class \code{teigen} or \code{mclust}.
}
\item{bestmodel}{
Results from the best model on the selected variable subset; an object of class \code{teigen} or \code{mclust}.
}
\item{chosenrelation}{
Numeric indication of the relationship chosen according to total model uncertainty. The number corresponds to exponent in the relationship: for instance, a value of '4' suggests the quartic relationship. If the value \code{"Full dataset"} is given, then the unreduced data provides the best model uncertainty; can be avoided by specifying \code{forcereduction=TRUE} in the function call.
}
\item{uncertainty}{
Total model uncertainty associated with the best relationship.
}
\item{allmodelfit}{
List containing the results (\code{teigen} or \code{mclust} objects) from the post-selection analysis on each variable subset. Number corresponds to the exponent in the relationship. For instance, \code{vscc_object$allmodelfit[[1]]} gives the results from the analysis on the variables selected by the linear relationship.
}
}
\references{
See \code{citation("vscc")} for the variable selection references. See also \code{citation("teigen")} and \code{citation("mclust")} if using those families of models via the \code{automate} call.
}
\author{
Jeffrey L. Andrews, Paul D. McNicholas
}


\seealso{
\code{\link{teigen}}, \code{\link{Mclust}}
}
\examples{
require("mclust")
data(banknote)
head(banknote)
bankrun <- vscc(banknote[,-1])
head(bankrun$topselected) #Show preview of selected variables
table(banknote[,1], bankrun$initialrun$classification) #Clustering results on full data set
table(banknote[,1], bankrun$bestmodel$classification) #Clustering results on reduced data set
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
