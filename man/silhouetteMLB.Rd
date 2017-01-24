\name{silhouetteMLB}
\alias{silhouetteMLB}
\title{
Silhouette Analysis
}
\description{
Analyze the clustering quality and generate a Silhouette Plot
}
\usage{
silhouetteMLB(data, fac, method = "cosine", plot = TRUE)
}
\arguments{
  \item{data}{ numeric matrix }
  \item{fac}{ clustering factor }
  \item{method}{ method to be used as distance function }
  \item{plot}{ logic, if TRUE plot a silhouette plot during the analysis }
}
\value{
numeric vector including the silhouette values of the data poointts in the input matrix
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
tmp <- getTestRunArgs(4)
silhouetteMLB(data = tmp$data, fac = tmp$fac)
}
