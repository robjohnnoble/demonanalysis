% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/speedestimates.R
\name{trans_rate}
\alias{trans_rate}
\title{Moran process transition rates}
\usage{
trans_rate(i, j, r, K)
}
\arguments{
\item{i}{initial population of the focus type}

\item{j}{new population of the focus type (should be i, i-1 or i+1)}

\item{r}{relative fitness of the focus type}

\item{K}{total population size}
}
\value{
The transition rate.
}
\description{
Moran process transition rates
}
\examples{
trans_rate(1, 2, 1.1, 4)
}
