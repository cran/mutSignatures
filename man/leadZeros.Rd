\name{leadZeros}
\alias{leadZeros}
\title{
Add Leading Zeros to a Number
}
\description{
Converts a numeric to character and adds leading zeros, in order to allow consistent sorting. Currently, this function is not used by the core framework functions anymore.
}
\usage{
leadZeros(n, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ number (integer) to be transformed }
  \item{m}{ maximum number that will be used to define how many leading zeros to attach }
}
\value{ input number "n" convverted to character and modified by attaching the required number of leading zeros}
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
leadZeros(15, 99999)
}
