% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demonanalysis.R
\name{plot_cum_dist}
\alias{plot_cum_dist}
\title{Plot cumulative density of variant allele frequencies versus inverse allele frequency}
\usage{
plot_cum_dist(file_or_dataframe, generation = NA, max_y = NA, ...)
}
\arguments{
\item{file_or_dataframe}{file or data frame containing columns "Frequency" and "Count"}

\item{generation}{Generation at which to make the measurement (default NA corresponds to the final Generation)}

\item{max_y}{max value of y-axis}

\item{...}{other parameters passed to plot}
}
\value{
plot displyed on screen
}
\description{
Plot cumulative density of variant allele frequencies versus inverse allele frequency
}
\examples{
plot_cum_dist(system.file("extdata", "output_allele_counts.dat", 
package = "demonanalysis", mustWork = TRUE))
}
