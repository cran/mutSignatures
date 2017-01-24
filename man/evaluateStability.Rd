\name{evaluateStability}
\alias{evaluateStability}
\title{Evaluate Results Stability}
\description{
Perform a final Stability check comparing the results from all iterations of the analysis
}
\usage{
evaluateStability(wall, hall, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wall}{numeric matrix including the w results from all the iterations of the analysis}
  \item{hall}{numeric matrix including the h results from all the iterations of the analysis}
  \item{params}{object (list) including all the parameters required for running tha analysis}
}
\details{
The function evvaluates the results from all iterations by performin a silhouette check. A silhouette plot will also be plotted.
}
\value{
list including all results from the stability checks. This includes the most likely signatures (centroids) and exposures. All information for plotting the silhoueette plot will also be returned}
\references{
More info at: \url{http://www.biotechworld.it/bioinf/2017/01/22/mutsignatures/}
}
\author{
Damiano Fantini
}
\note{
Affiliation: Meeks Lab, Department of Urology, Feinberg School of Medicine, Northwestern University, Chicago, IL 60611 USA
Software: Original MATLAB framework was described here: Cell Rep. 2013 Jan 31;3(1):246-59
}

\examples{
tmp <- getTestRunArgs(5)
evaluateStability(wall = tmp$W, hall = tmp$H, params = tmp$params)
}
