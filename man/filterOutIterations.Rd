\name{filterOutIterations}
\alias{filterOutIterations}
\title{
Remove Iterations that Generated Outlier Results
}
\description{
Core function called from within a deconvoluteMutCounts() call}
\usage{
filterOutIterations(wall, hall, cnt.errors, cnt.reconstructed, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wall}{ numeric matrix combining w results from all iterations }
  \item{hall}{ numeric matrix combining h results from all iterations }
  \item{cnt.errors}{ numeric matrix combining all residuals from all iterations }
  \item{cnt.reconstructed}{numeric matrix combining fitted values from all iterations}
  \item{params}{ object (list) including alll parameters for running the analysis }
}
\value{ list including all data required for running the subsequent stability check }
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
# No examples are provided for this function
message("no examples")
}
