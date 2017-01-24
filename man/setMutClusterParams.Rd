\name{setMutClusterParams}
\alias{setMutClusterParams}
\title{
Set Parameters for Extracting Mutational Signatures }
\description{
Create an object (list) including all parameters required for running the mutSignatures framework 
}
\usage{
setMutClusterParams(num.processes.toextract = 2, tot.iterations = 4, 
tot.cores = 5, remove.weak.muttypes = 0.01, 
remove.last.percent = 0.07, process.distance = "cosine", 
tot.Replicates = 100, eps = 2.2204e-16, 
stopconv = 10000, niter = 1e+06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{num.processes.toextract}{ number of signatures to extract }
  \item{tot.iterations}{ number of iterations per core}
  \item{tot.cores}{ number of cores to use for the analysis}
  \item{remove.weak.muttypes}{ threshold for filtering out under-represented mutation types }
  \item{remove.last.percent}{ threshold for removing outlier iteration results }
  \item{process.distance}{ method for calculating distances. Default method is 'cosine' }
  \item{tot.Replicates}{ number of replicates while checking stability }
  \item{eps}{ close-to-zero numeric value for replacing zeros and preventing negative values in the matrix }
  \item{stopconv}{ number of non-changing iterations before stopping NMF }
  \item{niter}{ total number of terations for NMF }
}
\value{
Object (list) including all parameters for running the analysis
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
setMutClusterParams()
}
