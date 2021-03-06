% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s03_all_functions.R
\name{extractSignatures}
\alias{extractSignatures}
\title{Extract Signatures from Genomic Mutational Catalogs.}
\usage{
extractSignatures(mutCountMatrix, params, bootStrap = TRUE)
}
\arguments{
\item{mutCountMatrix}{numeric matrix of mutation counts}

\item{params}{list including all parameters for performing the analysis}

\item{bootStrap}{logical, shall bootstrapping be performed}
}
\value{
list including the following elements
\enumerate{
   \item \bold{Wall}:  all extracted signatures
   \item \bold{Hall}: all extracted exposures
   \item \bold{mutCounts.reconstructed}: fitted values
   \item \bold{mutCounts.errors}: residuals
 }
}
\description{
Extract mutational signatures after the input Data and the 
input parameters have been checked andvalidated.
}
\details{
This is one of the core functions included in the original mutSignatures R library, 
and in the WTSI MATLAB framework. This is an internal function.
}
\examples{
x <- mutSignatures:::getTestRunArgs("extractSignatures")
y <- mutSignatures:::extractSignatures(mutCountMatrix = as.matrix(x$muts), 
                                       params = as.list(x$params), bootStrap = TRUE)
y$Wk[1:10,]


}
\references{
More information and examples about mutational signature analysis can be found here:
\enumerate{
  \item \bold{GitHub Repo}: \url{https://github.com/dami82/mutSignatures/}
  \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
  \item \bold{Sci Rep paper}, introducing mutS: \url{https://www.nature.com/articles/s41598-020-75062-0/}
  \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
  \item WTSI framework: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/}
 }
}
\author{
Damiano Fantini, \email{damiano.fantini@gmail.com}
}
\keyword{internal}
