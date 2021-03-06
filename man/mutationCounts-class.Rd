% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s01_all_classes.R
\docType{class}
\name{mutationCounts-class}
\alias{mutationCounts-class}
\alias{initialize,mutationCounts-method}
\title{Class mutationCounts.}
\usage{
\S4method{initialize}{mutationCounts}(.Object, x, muts, samples)
}
\arguments{
\item{.Object}{the mutationCounts object being built}

\item{x}{data.frame including mutation count values for each biological sample}

\item{muts}{data.frame including information about mutation types}

\item{samples}{data.frame including information about sample identifiers (unique names)}
}
\description{
Class mutationCounts defines objects storing Mutation COunts data.
}
\section{Slots}{

\describe{
\item{\code{counts}}{data.frame including information about mutation counts}

\item{\code{mutTypes}}{data.frame including information about mutation types}

\item{\code{sampleId}}{data.frame including information about sample identifiers}
}}

\author{
Damiano Fantini \email{damiano.fantini@gmail.com}
}
