\name{decipherMutationalProcesses}
\alias{decipherMutationalProcesses}
\title{Decipher Mutational Processes Contributing to a Collection of Genomic Mutations}
\description{Cancer cells accumulate DNA mutations as result of DNA damage and DNA repair processes. This computational framework allows to decipher mutational processes from cancer-derived somatic mutational catalogs }
\usage{
decipherMutationalProcesses(input, params)
}
\arguments{
  \item{input}{object (list) including a Mutation Count matrix. This is the output of a setMutCountObject() call}
  \item{params}{object (list) including the parameters required for running the analysis. This is the output of a setMutClusterParams() call}
}
\value{
list including all results of the analysis. The extracted signatures (processes) are included in the "processes" element of the list. The relative contribution of each signature in each sample is summarized in the "exposures" element}
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
decipherMutationalProcesses(input = tmp$mutCount.obj, params = tmp$params)
}
}
