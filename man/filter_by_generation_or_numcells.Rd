% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datacollation.R
\name{filter_by_generation_or_numcells}
\alias{filter_by_generation_or_numcells}
\title{Filter a data frame by Generation or NumCells}
\usage{
filter_by_generation_or_numcells(df, path, generation = NA,
  numcells = NA, num_parameters = NA)
}
\arguments{
\item{df}{data frame}

\item{path}{folder containing output.dat (needed for numcells only)}

\item{generation}{Generation at which to filter (default NA corresponds to no filtering)}

\item{numcells}{Number of cells at which to filter (default NA corresponds to no filtering)}

\item{num_parameters}{Number of parameters, accounting for the first set of columns in the dataframe; required if df represents multiple simulations}
}
\value{
the combined dataframe
}
\description{
Filter a data frame by Generation or NumCells
}
\details{
If both \code{generation} and \code{numcells} are provided then \code{numcells} 
takes precedent. If \code{numcells} is provided and \code{df} lacks a \code{NumCells} 
column then a \code{NumCells} column will be added (using the \code{output.dat} file 
in the folder specified by \code{path}), unless \code{df} contains multiple \code{seed} 
values, in which case an error will result.
}
\examples{
df <- read_delim_special(system.file("extdata", "output_allele_counts.dat", 
package = "demonanalysis", mustWork = TRUE))
filter_by_generation_or_numcells(df, NA, generation = 10, numcells = NA)
my_path <- system.file("extdata", package = "demonanalysis", mustWork = TRUE)
filter_by_generation_or_numcells(df, my_path, generation = NA, numcells = 100)
}
