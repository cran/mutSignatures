\name{do.nmf}
\alias{do.nmf}
\title{
Perform Non-negative Matrix Factorization
}
\description{
Perform Non-negative Matrix Factorization
}
\usage{
do.nmf(v, r, params, verbose = TRUE)
}
\arguments{
  \item{v}{
numeric matrix of Mutation Type Counts
}
  \item{r}{
numeric, number of signatures to extract
}
  \item{params}{
object (list) including all paramaters for running the analysis
}
  \item{verbose}{
logic, TRUE for printing information abouut the progress of the analysis
}
}
\value{
list including two elements. 
\item{w}{extracted signatures}
\item{h}{contribution of each signature in all the samples of the input mut count matrix}
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
tmp <- getTestRunArgs(3)
do.nmf(v = tmp$v, r = tmp$r, params = tmp$params)
}
