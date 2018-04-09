demonanalysis
========

Analysis and plotting of demon data

### Updating

To pull the latest github version into an Euler folder (using terminal):

```
cd /cluster/work/bewi/members/lebidm_nobelr/demonanalysis/
git pull https://github.com/robjohnnoble/demonanalysis
```

To install that version on an Euler profile (using R):

``` r
library(devtools)
install("/cluster/work/bewi/members/lebidm_nobelr/demonanalysis")
```

### Launch R on Euler

```
module load new
module load r/3.4.0
R
setwd("/cluster/work/bewi/members/lebidm_nobelr/demon")
```

### Set up

``` r
library(demonanalysis)

# set parameter values:
subfolder_name <- Mar20_Batch1
pars <- c("migration_edge_only", "mu_driver_birth", "seed") # parameters that were varied within the batch
final_values <- c(1, 1, 9) # maximum values of the parameters that were varied
start_size_range <- 500 + (0:5) * 1000 # NumCells at time of initial measurement for forecasting
gap_range <- (1:10)/10 # gap between time of initial measurement and second measurement
final_size <- 1E6 # waiting time is measured until tumour reaches this NumCells value

input_dir <- paste0("all_results/", subfolder_name) # folder containing results of a batch of simulations
output_dir_plots <- paste0("plots/", subfolder_name) # folder to receive image files
output_dir_data <- paste0("data/", subfolder_name) # folder to receive data files
```

### Create plots

``` r
create_plots_batch(input_dir, output_dir_plots, pars, final_values, type = "chart") # create plots
```

### Get complete data

``` r
data <- all_output(input_dir, pars, final_values) # combined data for a batch of simulations
data <- add_relative_time(data, start_size = 5500) # add columns useful for plotting trajectories
# check number of replicates per parameter set:
count_seeds(data)
```

### Get summary data

``` r
summary <- get_summary(data, start_size_range, gap_range, final_size) # summary data for each simulation, for each combination of gap and final_size
# check number of replicates per parameter set:
count_seeds(summary)
# check number of rows:
dim(summary)[1] == sum(count_seeds(summary)) * length(start_size_range) * length(gap_range)

cor_summary <- get_cor_summary(summary, c("DriverDiversity", "DriverEdgeDiversity")) # summary dataframe of correlations with "outcome"
wait_cor_summary <- get_wait_cor_summary(summary, c("DriverDiversity", "DriverEdgeDiversity")) # summary dataframe of correlations with "waiting_time"
```

### Write data

``` r
write.csv(data, paste0(output_dir_data, "data.csv"))
write.csv(summary, paste0(output_dir_data, "summary.csv"))
write.csv(cor_summary, paste0(output_dir_data, "cor_summary.csv"))
write.csv(wait_cor_summary, paste0(output_dir_data, "wait_cor_summary.csv"))
```

### Read data

``` r
data <- read.csv(paste0(output_dir_data, "data.csv"))
summary <- read.csv(paste0(output_dir_data, "summary.csv"))
cor_summary <- read.csv(paste0(output_dir_data, "cor_summary.csv"))
wait_cor_summary <- read.csv(paste0(output_dir_data, "wait_cor_summary.csv"))
```

### Copy plots from Euler

```
cd /Users/rnoble/Documents/MontpellierDocuments/Models/Demon/EulerPlots
scp -r rnoble@euler.ethz.ch:/cluster/work/bewi/members/lebidm_nobelr/demon/plots/* ./
```



