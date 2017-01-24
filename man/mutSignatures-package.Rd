\name{mutSignatures-package}
\alias{mutSignatures-package}
\alias{mutSignatures}
\docType{package}
\title{Computational Framework for Deciphering Cancer Mutational Signatures}
\description{Accumulation of DNA mutations is a hallmark of cancer. Cancer cells accumulate DNA mutations as result of DNA damage and DNA repair processes. This computational framework allows to characterize mutational signatures from cancer-derived somatic mutational catalogs}
\details{Please, cite the following article when using this software: Fantini D et al, Mol Biol Cell. 2017 Jan 1;28(1):192-200. 
More info at: http://www.biotechworld.it/bioinf/2017/01/22/mutsignatures/
Software developed in the Meeks Lab, Dept. of Urology, Northwestern University, Chicago, USA
This is a R implementation of the original MATLAB framework described here: Cell Rep. 2013 Jan 31;3(1):246-59.
}
\author{
Damiano Fantini, Joshua J Meeks
Maintainer: Damiano Fantini <damiano.fantini@gmail.com>
}
\examples{
  citation()
  \dontrun{
  tmp <- getTestRunArgs()
  result <- extractSignatures(mutCountMatrix = tmp$mutCount.obj$mutCounts, params = tmp$params)
  result
  plotSignatureExposures(results)
  }
}
