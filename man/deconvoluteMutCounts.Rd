\name{deconvoluteMutCounts}
\alias{deconvoluteMutCounts}
\title{Deconvolute Mutation Counts}
\description{Characterize mutational signatures from cancer-derived somatic mutational catalogs}
\usage{
deconvoluteMutCounts(input.mutCounts, params)
}
\arguments{
  \item{input.mutCounts}{
  numeric matrix of Mutation Type Counts
}
  \item{params}{
object (list) including all parameters required for running the analysis
}
}
\value{
list including all the results from the deconvolution analysis. This function is called within the decipherMutationalProcesses() function after parameters and input data have been validated.  
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
\dontrun{
tmp <- getTestRunArgs(2)
tmp
deconvoluteMutCounts(input.mutCounts = tmp$mutCount.obj$mutCounts[1:5,1:6], params = tmp$params)
}
}
