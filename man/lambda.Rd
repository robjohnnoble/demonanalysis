% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/speedestimates.R
\name{lambda}
\alias{lambda}
\title{Expected migration rate per deme}
\usage{
lambda(i, K, d, m, r1, r2, selection = "birth")
}
\arguments{
\item{i}{number of cells in the deme}

\item{K}{deme carrying capacity}

\item{d}{migration distance relative to 1/sqrt(K)}

\item{m}{migration rate per cell, relative to birth rate}

\item{r1}{reference birth rate}

\item{r2}{birth rate of migrating cells, relative to r1}

\item{selection}{if selection = "birth & migration" then the rate will be multiplied by a factor of r2 * r1}
}
\value{
A number, corresponding to the rate.
}
\description{
Expected migration rate per deme
}
\details{
Assumes that min(d / sqrt(K), 1) is the probability that a migration attempt will land outside the deme. 
Importantly, does not account for the chance that the migrating cell will land in an already-occupied deme.
}
\examples{
lambda(1, 2, 1, 0.1, 1, 1.1)
}