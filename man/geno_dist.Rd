% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sweepmetrics.R
\name{geno_dist}
\alias{geno_dist}
\title{Get genotype frequencies in increasing order of size}
\usage{
geno_dist(pop_df, generation = NA)
}
\arguments{
\item{pop_df}{Dataframe with column names "Identity", "Population" and "Generation"}

\item{generation}{Generation at which to get frequencies (default NA corresponds to final generation)}
}
\value{
vector of genotype frequencies
}
\description{
Get genotype frequencies in increasing order of size
}
\examples{
library(dplyr)
library(ggmuller)
phylo <- filter(driver_phylo, CellsPerSample == -1)
pop_df <- get_population_df(phylo)
dist <- geno_dist(pop_df)
barplot(dist, ylim = c(0, 1))
}
