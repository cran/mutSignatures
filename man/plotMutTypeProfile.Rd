% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s03_all_functions.R
\name{plotMutTypeProfile}
\alias{plotMutTypeProfile}
\title{Plot Mutation Signature Profiles.}
\usage{
plotMutTypeProfile(
  mutCounts,
  mutLabs,
  freq = TRUE,
  ylim = "auto",
  ylab = "Fraction of Variants",
  xlab = "Sequence Motifs",
  xaxis_cex = 0.475,
  cols = c("#4eb3d3", "#040404", "#b30000", "#bdbdbd", "#41ab5d", "#dd3497"),
  main = "MutType Profile"
)
}
\arguments{
\item{mutCounts}{data.frame including mutation types counts or frequencies, such as a
data.frame of mutation counts from samples, or mutation type frequencies from a mutational signature.}

\item{mutLabs}{character vector, labels to be used for the mutation types}

\item{freq}{logical, shall frequency be plotted rather than counts. Defaults to TRUE}

\item{ylim}{values used for ylim. Defaults to "auto" (ylim automatically set)}

\item{ylab}{string, used as y-axis title. Defaults to "Fraction of Variants"}

\item{xlab}{string, used as x-axis title. Defaults to "Sequence Motifs"}

\item{xaxis_cex}{numeric, cex value for the xaxis}

\item{cols}{character vector, indicates the colors to be used for the bars. It typically requires 6 colors.}

\item{main}{string, tutle of the plot. Defaults to "MutType Profile"}
}
\value{
NULL. A plot is printed to the active device.
}
\description{
Build a barplot to visualize the relative abundance of mutation counts in a mutational
signature or biological sample of interest.
}
\details{
This function is part of the user-interface set of tools included in mutSignatures. This is an exported function.
}
\references{
More information and examples about mutational signature analysis can be found here:
\enumerate{
  \item \bold{GitHub Repo}: \url{https://github.com/dami82/mutSignatures/}
  \item \bold{More info and examples} about the mutSignatures R library: \url{https://www.data-pulse.com/dev_site/mutsignatures/}
  \item \bold{Sci Rep paper}, introducing mutS: \url{https://www.nature.com/articles/s41598-020-75062-0/}
  \item \bold{Oncogene paper}, Mutational Signatures Operative in Bladder Cancer: \url{https://www.nature.com/articles/s41388-017-0099-6}
 }
}
\author{
Damiano Fantini, \email{damiano.fantini@gmail.com}
}
