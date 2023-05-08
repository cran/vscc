\name{vsccmanly}
\alias{vsccmanly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variable Selection for Skewed Clustering and Classification
}
\description{
Performs variable selection under a clustering framework. Accounts 
for mixtures of non-Gaussian distributions via the ManlyTransform (via 'ManlyMix'). 
}
\usage{
vsccmanly(x, G=2:9, numstart=100, selection="backward",forcereduction=FALSE,
                     initstart="k-means", seedval=2354)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data frame or matrix to perform variable selection on
}
 \item{G}{
Vector for the number of groups to consider during initialization and/or post-selection analysis. Default is 2-9. 
}
  \item{numstart}{
Number of random starts.
}
  \item{selection}{
Forward or backward transformation parameter selection. User may also choose to fit a full Manly mixture (options are 'forward', 'backward', or 'none').
}
  \item{forcereduction}{
Logical indicating if the full data set should be considered (FALSE) when selecting the `best' variable subset via total model uncertainty.
}
  \item{initstart}{
Method for initial starting values (options are 'k-means' or 'hierarchical').
}
  \item{seedval}{
Value of seed, used for k-means initialization.
}
}
\value{
\item{selected}{
A list containing the subsets of variables selected for each relation. Each set is numbered according to the number in the exponential of the relationship. For instance, \code{vscc_object$selected[[3]]} corresponds to the variable subset selected by the cubic relationship.
}
\item{wss}{
The within-group variance associated with each variable from the full data set.
}

\item{topselected}{
The best variable subset according to the total model uncertainty.
}
\item{initialrun}{
Results from the initial model, prior to variable selection; an object of class \code{ManlyMix}.
}
\item{bestmodel}{
Results from the best model on the selected variable subset; an object of class \code{ManlyMix}.
}
\item{variables}{
Variables used to fit the final model.
}
\item{chosenrelation}{
Numeric indication of the relationship chosen according to total model uncertainty. The number corresponds to exponent in the relationship: for instance, a value of '4' suggests the quartic relationship. If the value \code{"Full dataset"} is given, then the unreduced data provides the best model uncertainty; can be avoided by specifying \code{forcereduction=TRUE} in the function call.
}
\item{uncertainty}{
Total model uncertainty associated with the best relationship.
}
\item{allmodelfit}{
List containing the results (\code{ManlyMix} objects) from the post-selection analysis on each variable subset. Number corresponds to the exponent in the relationship. For instance, \code{vscc_object$allmodelfit[[1]]} gives the results from the analysis on the variables selected by the linear relationship.
}
}
\references{
See \code{citation("vscc")} for the variable selection references.
}
\author{
Jeffrey L. Andrews, Mackenzie R. Neal, Paul D. McNicholas
}

\seealso{
\code{\link{vscc}}
}
\examples{
\dontrun{
data(ais)
X=ais[,3:13]
aisfor=vsccmanly(as.data.frame(scale(X)),G=2:9,selection = "forward", forcereduction = TRUE,
                        initstart = "k-means",seedval=2354) 
aisfor$variables #Show selected variables
table(ais[,1], aisfor$bestmodel$id) #Clustering results on reduced data set
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
