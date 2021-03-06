% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demonanalysis.R
\name{create_plots_batch}
\alias{create_plots_batch}
\title{Create image files for every simulation in a batch}
\usage{
create_plots_batch(input_dir, type = "plot", file_type = "png",
  output_dir = NA, ...)
}
\arguments{
\item{input_dir}{base input directory name}

\item{type}{what type of images to create: "plot", "chart" or "origintimes" (or a vector containing two or more of these strings)}

\item{file_type}{either "pdf" or "png" (other values default to "pdf")}

\item{output_dir}{folder in which to save the image files}

\item{...}{additional arguments passed to the plotting function}
}
\value{
a set of image files
}
\description{
Create image files for every simulation in a batch
}
