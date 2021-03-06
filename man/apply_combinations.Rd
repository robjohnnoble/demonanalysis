% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demonanalysis.R
\name{apply_combinations}
\alias{apply_combinations}
\title{Apply a function to every combination of some sequences}
\usage{
apply_combinations(vec, fn, ...)
}
\arguments{
\item{vec}{vector of final values of the sequences (initial values are always zero)}

\item{fn}{function to apply to the values}

\item{...}{other arguments passed to fn}
}
\value{
result of applying fn to every combination of vec values
}
\description{
Apply a function to every combination of some sequences
}
\examples{
apply_combinations(c(2, 3), mean)
}
