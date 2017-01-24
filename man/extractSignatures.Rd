\name{extractSignatures}
\alias{extractSignatures}
\title{ Extract Signatures from Genomic Mutational Catalogs}
\description{
Extract mutational signatures after the input Data and the input parameters have been checked and validated
}
\usage{
extractSignatures(mutCountMatrix, params, verbose = TRUE)
}
\arguments{
  \item{mutCountMatrix}{numeric matrix of mutation counts}
  \item{params}{object (list) including all parameters for performing the analysis}
  \item{verbose}{logic, if TRUE the progress of the analysis will be printed to screen}
}
\value{
list including the following elements
\item{Wall}{extracted signatures} 
\item{Hall}{extracted effects} 
\item{mutCounts.reconstructed}{fitted values} 
\item{mutCounts.errors}{residuals} 
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
extractSignatures(mutCountMatrix = tmp$mutCount.obj$mutCounts[1:10,1:3], params = tmp$params)
}
