\name{bootstrapCancerGenomes}
\alias{bootstrapCancerGenomes}
\title{Bootstrap a Mutation Count Matrix}
\description{
Rearrange a Mutation count Matrix using the multivariate normal distribution. The function returns a bootstrapped Mutation Count matrix whose dimensions are identical to the input matrix
}
\usage{
bootstrapCancerGenomes(genomes)
}
\arguments{
  \item{genomes}{
a numeric matrix of Mutation Counts. Rows correspond to Mutation Types, columns to different samples
}
}
\value{
a numeric matrix of Mutation Counts. Rows correspond to Mutation Types, columns to different samples
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
