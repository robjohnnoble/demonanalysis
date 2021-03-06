% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demonanalysis.R
\name{hist2}
\alias{hist2}
\title{Alternative to base hist function (using dplyr)}
\usage{
hist2(x, breaks, counts = 1)
}
\arguments{
\item{x}{a vector of values for which the histogram is desired}

\item{breaks}{a vector giving the breakpoints between histogram bins}

\item{counts}{optional vector of counts for each x value}
}
\value{
data frame with counts and densities
}
\description{
Alternative to base hist function (using dplyr)
}
\examples{
freq <- seq(0, 1, length = 20)
num <- rbinom(20, 10, 0.5)
breaks <- seq(0, 1, length.out = 5)
hist2(freq, breaks, num)
hist2(freq, breaks)

# equivalent using standard function:
hist(rep(x = freq, times = num), breaks, plot = FALSE)
hist(freq, breaks, plot = FALSE)
}
