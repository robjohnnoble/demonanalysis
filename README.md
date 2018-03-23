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

### Basic usage

``` r
library(demonanalysis)

input_dir <- "all_results/Mar20_Batch1" # folder containing results of a batch of simulations
output_dir <- "plots/Mar20_Batch1" # folder to receive image files
pars <- c("migration_edge_only", "mu_driver_birth", "seed") # parameters that were varied within the batch
final_values <- c(1, 1, 9) # maximum values of the parameters that were varied

create_plots_batch(input_dir, output_dir, pars, final_values) # create the plots
```
