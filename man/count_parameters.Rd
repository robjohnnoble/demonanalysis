% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datacollation.R
\name{count_parameters}
\alias{count_parameters}
\title{Count the number of parameters in a parameter file in a specified folder or its subfolders.}
\usage{
count_parameters(full_dir)
}
\arguments{
\item{full_dir}{folder name}
}
\value{
number of parameters
}
\description{
Count the number of parameters in a parameter file in a specified folder or its subfolders.
}
\examples{
count_parameters(system.file("extdata", "", package = "demonanalysis", mustWork = TRUE))
count_parameters(system.file("example_batch", "", package = "demonanalysis", mustWork = TRUE))
}
