\name{addWeak}
\alias{addWeak}
\title{Add Weak Mutation Types}
\description{Restore Mutation Types that were initially excluded because a low number of total counts}
\usage{
addWeak(mutationTypesToAddSet, processes_I, processesStd_I, 
Wall_I, genomeErrors_I, genomesReconstructed_I)
}
\arguments{
  \item{mutationTypesToAddSet}{
Set of mutations to restore
}
  \item{processes_I}{
Set of Mutational Processes
}
  \item{processesStd_I}{
Set of standard deviations of all Mutational Processes
}
  \item{Wall_I}{
Set of all W matrices previously extracted
}
  \item{genomeErrors_I}{
Set of all residuals
}
  \item{genomesReconstructed_I}{
Fitted Values according to the most likely Model
}
}
\value{
Output is the final result of the deconvolution process
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
# there are no examples for this function
message("No examples")
}