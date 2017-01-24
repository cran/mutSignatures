\name{getTestRunArgs}
\alias{getTestRunArgs}
\title{Generate Arguments for Running Examples and Mock Runs}
\description{
this function generates objects that can be used for running the examples in the package documentation files as well as some simple runs
}
\usage{
getTestRunArgs(testN = 1)
}
\arguments{
  \item{testN}{ numeric, index corresponding to the type of test input to generate }
}
\value{
list of objects to be used as arguments for running the example codes of the package }
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
getTestRunArgs(1)
}
