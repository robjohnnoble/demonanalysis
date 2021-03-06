% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datacollation.R
\name{corCI_low}
\alias{corCI_low}
\title{Function to return the lower bound of the confidence interval of a Spearman's rank correlation coefficient.
This function is created to work with find_correlations}
\usage{
corCI_low(Var1, Var2)
}
\arguments{
\item{Var1}{first variable of interest}

\item{Var2}{second variable of interest}
}
\value{
lower bound of the confidence interval of a Spearman's rank correlation coefficient, computed by bootpstraping, 
with nrep = 500 the number of replicates for bootstraping and conf.level = 0.95 the confidence level of the interval.
return NA if no confidence interval has been computed, instead of NULL.
}
\description{
Function to return the lower bound of the confidence interval of a Spearman's rank correlation coefficient.
This function is created to work with find_correlations
}
