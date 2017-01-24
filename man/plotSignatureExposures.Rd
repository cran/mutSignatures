\name{plotSignatureExposures}
\alias{plotSignatureExposures}
\title{Generate Plot Signature}
\description{Plot Signature Contributions of All Samples in a Dataset }
\usage{
plotSignatureExposures(mutSign.analysis.result)
}
\arguments{
  \item{mutSign.analysis.result}{ object (list) including the results of a decipherMutationalProcesses() analysis }
}
\details{Data are plotted using ggplot2}
\value{ ggplot2 object }
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
tmp <- getTestRunArgs(6)
gp <- plotSignatureExposures(tmp)
gp + ggtitle("de novo mutSignatures")  + scale_fill_manual(values = c("#377eb8", "#4daf4a"))
}
