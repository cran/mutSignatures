\name{setMutCountObject}
\alias{setMutCountObject}
\title{
Create a Mutation Count Object Suitable for Signatures Extraction
}
\description{
Create a Mutation Count Object that can be used as data input for running the Signatures Extraction Framework
}
\usage{
setMutCountObject(mutCountMatrix, mutationTypes = NULL, 
sampleNames = NULL, datasetName = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutCountMatrix}{ numeric matrix of mutation counts }
  \item{mutationTypes}{ vector including the mutation types }
  \item{sampleNames}{ vector including sample names }
  \item{datasetName}{ name of the dataset }
}
\value{
list including the mutation Count Matrix and input data for de-novo extracting Mutational Signatures
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
tmp <- getTestRunArgs()[["mutCount.obj"]][["mutCounts"]]
tmp.obj <- setMutCountObject(tmp)
tmp.obj
}
