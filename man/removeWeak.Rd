\name{removeWeak}
\alias{removeWeak}
\title{Remove Mutation Types Not Meeting the Threshold }
\description{
Remove mutation types that account for a total number of mutations below a defined threshold.
}
\usage{
removeWeak(input.mutCounts, params)
}
\arguments{
  \item{input.mutCounts}{numeric matrix of Mutation Counts}
  \item{params}{object (list) including all parameters required for running the analysis}
}
\value{
list including two elements
  \item{removed.mutset}{numeric vector including the ids of the removed rows. If removed.mutset == -1, no rows were removed }
  \item{output.mutCounts}{filtered Mutation Count matrix  }
}
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
tmp <- getTestRunArgs()
removeWeak(input.mutCounts = tmp$mutCount.obj$mutCounts, params = tmp$params)
}
