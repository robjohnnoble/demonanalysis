demonanalysis
========

Analysis and plotting of demon data

## From a single simulation

To plot charts of variant allele frequencies and genotype sizes use `plot_all_charts`. For normal tissue, you should specify a value for `max_allele_count` so that the axes are appropriately scaled.

To plot Muller plots and grids use `plot_all_images`.

## From a batch of simulations run on Euler

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

subfolder_name <- "XXX" # insert name of top-level folder

input_dir <- paste0("all_results/", subfolder_name) # folder containing results of the batch
num_parameters <- count_parameters(input_dir) # number of simulation parameters (first columns in data)
output_dir_plots <- paste0("plots/", subfolder_name) # folder to receive image files
output_dir_data <- paste0("data/", subfolder_name) # folder containing data files
```

### Check which simulations completed successfully

``` r
all_statuses(input_dir, summary = TRUE) # should be "Exit code 0" (when finished) or "So far no status" (while running)
```

### Create plots for a batch of simulations

For growing tumours, to plot charts of variant allele frequencies and genotype sizes:

``` r
create_plots_batch(input_dir, output_dir = output_dir_plots, type = "chart")
```

For normal tissue, to plot charts of variant allele frequencies and genotype sizes:

``` r
create_plots_batch(input_dir, output_dir = output_dir_plots, type = "chart", max_genotype_size = 50, max_allele_count = 50)
```

To plot Muller plots and grids:

``` r
create_plots_batch(input_dir, output_dir = output_dir_plots, type = "plot")
```

### Get general-purpose data

``` r
data <- all_output(input_dir, include_diversities = FALSE) # combined data for a batch of simulations, excluding diversity columns
```

### Plot relationships between variables (using dataframe)

``` r
plot_curves_faceted(data, num_parameters, x_var = "Generation", y_var = "MeanBirthRate", output_filename = "curves", output_dir = output_dir_plots)
# change x_var, y_var and output_filename as appropriate
```

### Get additional data for forecasting

``` r
data <- all_output(input_dir) # combined data for a batch of simulations, including diversity columns
data <- add_relative_time(data, start_size = 5500, num_parameters = num_parameters) # add columns useful for plotting trajectories

start_size_range <- 500 * 2^(0:8) # NumCells at time of initial measurement for forecasting
gap_range <- (1:10)/10 # gap between time of initial measurement and second measurement
final_size <- 1E5 # waiting time is measured until tumour reaches this NumCells value

summary <- get_summary(data, start_size_range, gap_range, final_size, num_parameters = num_parameters) # summary data for each simulation, for each combination of gap and final_size

cor_summary <- get_cor_summary(summary, c("DriverDiversity", "DriverEdgeDiversity"), num_parameters = num_parameters, min_count = 5) # summary dataframe of correlations with "outcome", including all cells
wait_cor_summary <- get_wait_cor_summary(summary, c("DriverDiversity", "DriverEdgeDiversity"), num_parameters = num_parameters, min_count = 5) # summary dataframe of correlations with "waiting_time", including all cells
depth_wait_cor_summary <- get_wait_cor_summary(summary, c(paste0("DriverDiversityFrom1SamplesAtDepth", 0:10), paste0("DriverDiversityFrom4SamplesAtDepth", 0:10)), num_parameters, min_count = 5) # summary dataframe of correlations with "waiting_time" for different biopsy protocols
```

### Write data

``` r
write.csv(data, paste0(output_dir_data, "/data.csv"), row.names = FALSE)
```

### Read data

``` r
library(readr)
data <- read_csv(paste0(output_dir_data, "/data.csv"), guess_max = 1E4) # large value of guess_max improves guessing of column types
```

### Copy plots from Euler

```
cd /Users/rnoble/Documents/MontpellierDocuments/Models/Demon/EulerPlots/April_6th_batch1/Charts
scp -r rnoble@euler.ethz.ch:/cluster/work/bewi/members/lebidm_nobelr/demon/plots/April_6th_batch1/chart* ./
```



