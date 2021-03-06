% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s02_all_methods.R
\name{coerceObj}
\alias{coerceObj}
\alias{coerceObj,mutFrameworkParams,character-method}
\alias{coerceObj,mutationSignatures,character-method}
\alias{coerceObj,mutationCounts,character-method}
\alias{coerceObj,mutSignExposures,character-method}
\alias{coerceObj,data.frame,character-method}
\title{Method coerceObj.}
\usage{
coerceObj(x, to, ...)

\S4method{coerceObj}{mutFrameworkParams,character}(x, to)

\S4method{coerceObj}{mutationSignatures,character}(x, to)

\S4method{coerceObj}{mutationCounts,character}(x, to, ...)

\S4method{coerceObj}{mutSignExposures,character}(x, to, ...)

\S4method{coerceObj}{data.frame,character}(x, to, ...)
}
\arguments{
\item{x}{an object to coerce to a different format}

\item{to}{string, indicates the expected format (such as list or data.frame)}

\item{...}{additional parameters passed to the functions used for the coercion}
}
\description{
Cast an object to a different format, by extracting and returning the
most appropriate information. Note that data.frames can be coerced to one
of the classes defined in the mutSignatures package using coerceObj.
}
