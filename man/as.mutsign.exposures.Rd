% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s02_all_methods.R
\name{as.mutsign.exposures}
\alias{as.mutsign.exposures}
\alias{as.mutsign.exposures,data.frame,logical-method}
\title{Method as.mutsign.exposures.}
\usage{
as.mutsign.exposures(x, samplesAsCols = TRUE)

\S4method{as.mutsign.exposures}{data.frame,logical}(x, samplesAsCols = TRUE)
}
\arguments{
\item{x}{a data.frame to be converted to a mutSignExposures-class object.}

\item{samplesAsCols}{logical, are samples listed as columns in the input data.frame. If FALSE, 
samples are expected to be listed as rows in the input data.frame}
}
\description{
Cast a data.frame into a mutSignExposures-class object.
}
