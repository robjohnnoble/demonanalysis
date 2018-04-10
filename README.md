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

subfolder_name <- "April_6th_batch1" # batch name
start_size_range <- 500 + (0:5) * 1000 # NumCells at time of initial measurement for forecasting
gap_range <- (1:10)/10 # gap between time of initial measurement and second measurement
final_size <- 1E5 # waiting time is measured until tumour reaches this NumCells value

input_dir <- paste0("all_results/", subfolder_name) # folder containing results of the batch
num_parameters <- count_parameters(input_dir) # number of simulation parameters (first columns in data)
output_dir_plots <- paste0("plots/", subfolder_name) # folder to receive image files
output_dir_data <- paste0("data/", subfolder_name) # folder containing data files
```

### Create plots

``` r
create_plots_batch(input_dir, output_dir_plots, pars, final_values, type = "chart") # change type to "chart", "plot" or c("chart", "plot") as needed
```

### Get complete data

``` r
data <- all_output(input_dir) # combined data for a batch of simulations
data <- add_relative_time(data, start_size = 5500, num_parameters = num_parameters) # add columns useful for plotting trajectories
```

### Get summary data

``` r
summary <- get_summary(data, start_size_range, gap_range, final_size, num_parameters = num_parameters) # summary data for each simulation, for each combination of gap and final_size

cor_summary <- get_cor_summary(summary, c("DriverDiversity", "DriverEdgeDiversity"), num_parameters = num_parameters, min_count = 5) # summary dataframe of correlations with "outcome"
wait_cor_summary <- get_wait_cor_summary(summary, c("DriverDiversity", "DriverEdgeDiversity"), num_parameters = num_parameters, min_count = 5) # summary dataframe of correlations with "waiting_time"
```

### Write data

``` r
write.csv(data, paste0(output_dir_data, "/data.csv"), row.names = FALSE)
write.csv(summary, paste0(output_dir_data, "/summary.csv"), row.names = FALSE)
write.csv(cor_summary, paste0(output_dir_data, "/cor_summary.csv"), row.names = FALSE)
write.csv(wait_cor_summary, paste0(output_dir_data, "/wait_cor_summary.csv"), row.names = FALSE)
```

### Read data

``` r
data <- read.csv(paste0(output_dir_data, "/data.csv"))
summary <- read.csv(paste0(output_dir_data, "/summary.csv"))
cor_summary <- read.csv(paste0(output_dir_data, "/cor_summary.csv"))
wait_cor_summary <- read.csv(paste0(output_dir_data, "/wait_cor_summary.csv"))
```

### Copy plots from Euler

```
cd /Users/rnoble/Documents/MontpellierDocuments/Models/Demon/EulerPlots
scp -r rnoble@euler.ethz.ch:/cluster/work/bewi/members/lebidm_nobelr/demon/plots/* ./
```



