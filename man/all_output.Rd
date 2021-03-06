% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datacollation.R
\name{all_output}
\alias{all_output}
\title{Create a composite dataframe for a batch of simulations, derived from multiple 
data files per simulation.}
\usage{
all_output(input_dir, include_diversities = TRUE, df_type = "output",
  max_generation = FALSE, vaf_cut_off = NA, generation = NA,
  numcells = NA, n_cores = NA, ExitCode4 = FALSE)
}
\arguments{
\item{input_dir}{base input directory name}

\item{include_diversities}{boolean whether to include diversity metrics}

\item{df_type}{which dataframes to combine}

\item{max_generation}{If TRUE then results are returned only for the final generation}

\item{vaf_cut_off}{exclude genotypes with vaf lower cut off from combined_df}

\item{generation}{Generation at which to filter (default NA corresponds to no filtering)}

\item{numcells}{Number of cells at which to filter (default NA corresponds to no filtering)}

\item{n_cores}{Number of cores to use (default NA is not parallelized)}

\item{ExitCode4}{: if TRUE, this means that we want to include in the analysis the simulations whose error message is Exit Code 4. This 
will call the function combine_dfs with argument ExitCode4=TRUE.}
}
\value{
a combined dataframe
}
\description{
Create a composite dataframe for a batch of simulations, derived from multiple 
data files per simulation.
}
\details{
If both \code{generation} and \code{numcells} are provided then \code{numcells} 
takes precedent.
}
\examples{
all_output(system.file("example_batch", "", package = "demonanalysis", mustWork = TRUE))
all_output(system.file("example_batch", "", package = "demonanalysis", mustWork = TRUE), 
df_type = "driver_genotype_properties")
all_output(system.file("example_batch", "", package = "demonanalysis", mustWork = TRUE), 
df_type = "genotype_properties", vaf_cut_off = 0.002)
all_output(system.file("example_batch", "", package = "demonanalysis", mustWork = TRUE), 
df_type = "allele_counts", generation = 10)
}
