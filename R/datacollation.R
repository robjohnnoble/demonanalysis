#' Count the number of parameters in a parameter file in a specified folder or its subfolders.
#' 
#' @param full_dir folder name
#' 
#' @return number of parameters
#' 
#' @importFrom readr read_delim
#' 
#' @export
#' 
#' @examples
#' count_parameters(system.file("extdata", "", package = "demonanalysis", mustWork = TRUE))
#' count_parameters(system.file("example_batch", "", package = "demonanalysis", mustWork = TRUE))
count_parameters <- function(full_dir) {
  
  if(substr(full_dir, nchar(full_dir), nchar(full_dir)) == "/") full_dir <- substr(full_dir, 1, nchar(full_dir) - 1)
  
  while(1) {
    contents_list = list.files(full_dir)
    
    exists_list <- dir.exists(paste0(full_dir, "/", contents_list))
    dirs_list <- contents_list[exists_list]
    
    if("parameters.dat" %in% contents_list) {
      break
    } else if(length(dirs_list) > 0) {
      sub_dir <- dirs_list[1]
    } else {
      stop(paste0("Cannot find a parameters file in ", full_dir))
    }
    
    full_dir <- paste0(full_dir, "/", sub_dir)
  }
  
  df_pars <- read_delim(paste0(full_dir, "/parameters.dat"), "\t")
  
  return(dim(df_pars)[2])
}

#' Find the names and final values of parameters that were varied in a batch of simulations
#' 
#' @param input_dir base input directory name
#' 
#' @return dataframe of parameter names and values, or NA if input_dir contains results of a single simulation
#' 
#' @export
#' 
#' @examples
#' parameter_names_and_values(system.file("example_batch", "", 
#' package = "demonanalysis", mustWork = TRUE))
parameter_names_and_values <- function(input_dir) {
  
  if(substr(input_dir, nchar(input_dir), nchar(input_dir)) == "/") input_dir <- substr(input_dir, 1, nchar(input_dir) - 1)
  
  parent_dir <- input_dir
  
  if("output.dat" %in% list.files(parent_dir, recursive = FALSE, full.names = FALSE)) return(NA)
  
  out_df <- data.frame("name" = NULL, "final_value" = NULL)
  
  repeat{
    dirs_list <- list.dirs(parent_dir, recursive = FALSE, full.names = FALSE) # list of subfolders
    
    if(identical(dirs_list, character(0))) stop(paste0("Invalid folder name (", parent_dir), ")")
    
    list_splits <- strsplit(dirs_list, "_")
    values <- lapply(list_splits, function (x) x[length(x)])
    values <- as.numeric(unlist(values))
    final_dir <- dirs_list[which.max(values)] # final subfolder (with largest parameter value)
    
    splits <- strsplit(final_dir, "_")[[1]]
    parameter_val <- as.numeric(splits[length(splits)])
    parameter_name <- substr(final_dir, 1, nchar(final_dir) - nchar(parameter_val) - 1)
    
    out_df <- rbind(out_df, data.frame("name" = parameter_name, "final_value" = parameter_val))
    
    parent_dir <- paste0(parent_dir, "/", final_dir)
    
    if("parameters.dat" %in% list.files(parent_dir, recursive = FALSE, full.names = FALSE)) break
  }
  return(out_df)
}

#' Add derived variables to a dataframe
#' 
#' @param df dataframe with columns including "Generation" and "NumCells"
#' @param num_parameters number of parameters, accounting for the first set of columns in the dataframe
#' 
#' @return The same dataframe with additional columns. For each simulation, 
#' "maxgen" is the maximum value of Generation; "gen_adj" is the time elapsed
#' since Generation zero, relative to maxgen; SmoothNumCells is NumCells with loess smoothing; 
#' SmoothRadius is the radius of a disc of area NumCells, after loess smoothing; 
#' GrowthRate and RadiusGrowthRate are the rates of change of SmoothNumCells and SmoothRadius, 
#' relative to Generation.
#' 
#' @importFrom stats loess
#' 
#' @export
#' 
#' @examples
#' df <- data.frame(p = 1, seed = rep(1:2, each = 10), 
#' Generation = c(1:10, 3:12), NumCells = rep(1:10, times = 2) + rnorm(20, 0, 0.1))
#' add_columns(df, 2)
add_columns <- function(df, num_parameters) {
  df <- df %>% group_by_at(1:num_parameters) %>% 
    mutate(maxgen = max(Generation, na.rm = TRUE)) %>% 
    mutate(gen_adj = Generation / maxgen) %>% 
    ungroup()
  
  df <- df %>% group_by_at(1:num_parameters) %>% 
    mutate(SmoothNumCells = 10^loess(log10(NumCells) ~ log10(Generation + 1), span = 0.75)$fitted) %>% 
    mutate(SmoothRadius = 10^loess(log10(sqrt(NumCells/pi)) ~ log10(Generation + 1), span = 0.75)$fitted) %>% 
    ungroup()
  
  df <- df %>% group_by_at(1:num_parameters) %>% 
    mutate(GrowthRate = (SmoothNumCells - lag(SmoothNumCells, 1)) / (Generation - lag(Generation, 1))) %>% 
    mutate(RadiusGrowthRate = (SmoothRadius - lag(SmoothRadius, 1)) / (Generation - lag(Generation, 1))) %>% 
    ungroup()
  
  # replace NA values:
  df[is.na(df$GrowthRate), "GrowthRate"] <- df[!is.na(df$GrowthRate) & is.na(lag(df$GrowthRate, 1)), "GrowthRate"]
  df[is.na(df$RadiusGrowthRate), "RadiusGrowthRate"] <- df[!is.na(df$RadiusGrowthRate) & is.na(lag(df$RadiusGrowthRate, 1)), "RadiusGrowthRate"]

  return(df)
}

#' Add relative time column to a dataframe, and diversity at the new time zero
#' 
#' @param df dataframe with columns including "gen_adj", "NumCells" and "DriverEdgeDiversity"
#' @param start_size value of NumCells at which new_time should be set to zero
#' @param num_parameters number of parameters, accounting for the first set of columns in the dataframe
#' 
#' @return the same dataframe with additional columns: "new_time" is the time elapsed since
#' NumCells = start_size; "div0" is DriverEdgeDiversity when NumCells = start_size; 
#' "rank_div0" is div0 rescaled to a rank between 0 and 1
#' 
#' @importFrom scales rescale
#' 
#' @export
#' 
#' @examples
#' num_parameters = count_parameters(system.file("extdata", "", 
#' package = "demonanalysis", mustWork = TRUE))
#' comb_df <- combine_dfs(system.file("extdata", "", package = "demonanalysis", mustWork = TRUE))
#' add_relative_time(comb_df, start_size = 100, num_parameters = num_parameters)
add_relative_time <- function(df, start_size, num_parameters) {
  col_nums <- c(1:num_parameters)
  
  df <- df %>% group_by_at(col_nums) %>% 
    mutate(new_time = gen_adj - min(gen_adj[NumCells >= start_size], na.rm = TRUE)) %>% 
    mutate(div0 = min(DriverEdgeDiversity[gen_adj == min(gen_adj[NumCells >= start_size & 
           (!is.na(DriverEdgeDiversity) | Generation == max(Generation))], na.rm = TRUE)])) %>% 
    mutate(GrowthRate0 = min(GrowthRate[gen_adj == min(gen_adj[NumCells >= start_size & 
           (!is.na(GrowthRate) | Generation == max(Generation))], na.rm = TRUE)])) %>% 
    ungroup()
  
  col_nums <- col_nums[col_nums != which(colnames(df) == "seed")]
    
  df <- df %>% group_by_at(col_nums) %>% 
    mutate(rank_div0 = rescale(rank(div0))) %>% 
    mutate(rank_GrowthRate0 = rescale(rank(GrowthRate0))) %>% 
    ungroup()

  return(df)
}

#' Filter a data frame by Generation or NumCells
#' 
#' @param df data frame
#' @param path folder containing output.dat (needed for numcells only)
#' @param generation Generation at which to filter (default NA corresponds to no filtering)
#' @param numcells Number of cells at which to filter (default NA corresponds to no filtering)
#' @param num_parameters Number of parameters, accounting for the first set of columns in the dataframe; required if df represents multiple simulations
#' 
#' @return the combined dataframe
#' 
#' @export
#' 
#' @details If both \code{generation} and \code{numcells} are provided then \code{numcells} 
#' takes precedent. If \code{numcells} is provided and \code{df} lacks a \code{NumCells} 
#' column then a \code{NumCells} column will be added (using the \code{output.dat} file 
#' in the folder specified by \code{path}), unless \code{df} contains multiple \code{seed} 
#' values, in which case an error will result.
#' 
#' @examples
#' df <- read_delim_special(system.file("extdata", "output_allele_counts.dat", 
#' package = "demonanalysis", mustWork = TRUE))
#' filter_by_generation_or_numcells(df, NA, generation = 10, numcells = NA)
#' my_path <- system.file("extdata", package = "demonanalysis", mustWork = TRUE)
#' filter_by_generation_or_numcells(df, my_path, generation = NA, numcells = 100)
filter_by_generation_or_numcells <- function(df, path, generation = NA, numcells = NA, num_parameters = NA) {
  # count simulations:
  if(!is.na(num_parameters)) {
    df <- group_by_at(df, 1:num_parameters)
    num_sims <- n_groups(df)
    df <- ungroup(df)
  }
  else num_sims <- 1
  if(num_sims > 1 || ("seed" %in% colnames(df) && length(unique(df$seed)) > 1)) multiple_sims <- 1
  else multiple_sims <- 0
  
  if(!is.na(numcells)) {
    # add NumCells column if needed (when possible):
    if(!("NumCells" %in% colnames(df))) {
      if(multiple_sims) stop("Cannot add NumCells column to a dataframe that represents multiple simulations.")
      # add NumCells column from output.dat:
      if(substr(path, nchar(path), nchar(path)) == "/") path <- substr(path, 1, nchar(path) - 1)
      ref_df <- read_delim_special(paste0(path, "/output.dat"))
      ref_df <- select(ref_df, Generation, NumCells) %>% 
        filter(Generation %in% df$Generation)
      df <- merge(df, ref_df)
    }
    # filter by closest NumCells to user input:
    if(multiple_sims) {
      if(is.na(num_parameters)) stop("Need to specify num_parameters for data frame representing multiple simulations.")
      df <- group_by_at(df, 1:num_parameters) %>% 
        filter(abs(NumCells - numcells) == min(abs(NumCells - numcells))) %>% 
        filter(NumCells == min(NumCells)) %>% 
        ungroup()
    }
    else {
      print(paste0("About to filter; numcells = ", numcells))
      print("NumCells:")
      print(unique(df$NumCells))
      df <- filter(df, abs(NumCells - numcells) == min(abs(NumCells - numcells))) %>% 
        filter(NumCells == min(NumCells))
      print("Filtered")
    }
  }
  else if(!is.na(generation)) {
    # filter by closest Generation to user input:
    if(multiple_sims) {
      if(is.na(num_parameters)) stop("Need to specify num_parameters for data frame representing multiple simulations.")
      df <- group_by_at(df, 1:num_parameters) %>% 
        filter(abs(Generation - generation) == min(abs(Generation - generation))) %>% 
        filter(Generation == min(Generation)) %>% 
        ungroup()
    }
    else {
      df <- filter(df, abs(Generation - generation) == min(abs(Generation - generation))) %>% 
        filter(Generation == min(Generation))
    }
  }
  return(df)
}

#' Combine (and optionally process) data for one simulation
#' 
#' @param full_dir base input directory name
#' @param include_diversities boolean whether to include diversity metrics (if df_type == "output")
#' @param df_type which dataframes to combine
#' @param vaf_cut_off exclude genotypes with vaf lower than cut off (if df_type == "genotype_properties" or "driver_genotype_properties")
#' @param generation Generation at which to filter (default NA corresponds to no filtering)
#' @param numcells Number of cells at which to filter (default NA corresponds to no filtering)
#' @param num_parameters Number of parameters, accounting for the first set of columns in the dataframe (optional, but may improve speed)
#' 
#' @return the combined dataframe
#' 
#' @importFrom readr read_delim
#' @importFrom dplyr bind_cols
#' @importFrom data.table fread
#' @importFrom moments skewness
#' 
#' @export
#' 
#' @details If both \code{generation} and \code{numcells} are provided then \code{numcells} 
#' takes precedent.
#' 
#' @examples
#' combine_dfs(system.file("extdata", "", package = "demonanalysis", mustWork = TRUE))
#' combine_dfs(full_dir = system.file("extdata", "", package = "demonanalysis", mustWork = TRUE), 
#' df_type = "driver_genotype_properties")
#' combine_dfs(full_dir = system.file("extdata", "", package = "demonanalysis", mustWork = TRUE), 
#' df_type = "genotype_properties", vaf_cut_off = 0.002)
#' combine_dfs(full_dir = system.file("extdata", "", package = "demonanalysis", mustWork = TRUE), 
#' df_type = "genotype_counts", numcells = 100)
combine_dfs <- function(full_dir, include_diversities = TRUE, df_type = "output", vaf_cut_off = NA, generation = NA, numcells = NA, num_parameters = NA) {
  
  if(substr(full_dir, nchar(full_dir), nchar(full_dir)) == "/") full_dir <- substr(full_dir, 1, nchar(full_dir) - 1)
  
  file_pars <- paste0(full_dir, "/parameters.dat")
  file_out <- paste0(full_dir, "/output.dat")
  file_div <- paste0(full_dir, "/output_diversities.dat")
  file_driver_phylo <- paste0(full_dir, "/driver_phylo.dat")
  file_allele_counts <- paste0(full_dir, "/output_allele_counts.dat")
  
  df_out <- fread(file_out)
  df_pars <- fread(file_pars)
  
  if (df_type == "output"){
    # procedure for 'traditional' df_out (output.dat)
    if(include_diversities) df_div <- fread(file_div)
    df_driver_phylo <- fread(file_driver_phylo)
    
    df_driver_phylo <- filter(df_driver_phylo, CellsPerSample == -1, NumSamples == 1, SampleDepth == -1)
    df_driver_phylo <- df_driver_phylo[!duplicated(df_driver_phylo), ]
    pop_df <- get_population_df(df_driver_phylo)
    
    if(include_diversities) temp <- merge(df_out, df_div, all = TRUE)
    else temp <- df_out
    
    # add parameter columns:
    if(nrow(temp) == 0){
      temp <- NULL
    } else {
      temp <- cbind(df_pars, temp)
    }
    
    # adds maxgen and gen_adj columns
    if(is.na(num_parameters)) num_parameters <- count_parameters(full_dir)
    temp <- add_columns(temp, num_parameters) 
    
    # add sweep_seq columns (specific for output.dat?)
    sweep_seq <- sweep_sequence(pop_df, lag_type = "proportions", breaks = 10)
    temp <- mutate(temp, mean_autocor = mean(sweep_seq), 
                     log_mean_autocor = log(mean(sweep_seq)), 
                     sqrt_mean_autocor = sqrt(mean(sweep_seq)), 
                     skewness = skewness(sweep_seq))
  } else if (df_type %in% c("allele_counts", "driver_allele_counts", "genotype_counts", "driver_genotype_counts")){
    temp <- fread(paste0(full_dir, "/output_", df_type, ".dat"))
    # add parameter columns:
    if(nrow(temp) == 0){
      temp <- NULL
    } else {
      temp <- cbind(df_pars, temp)
    }
  } else if (df_type %in% c("driver_genotype_properties", "genotype_properties")){
    temp <- fread(paste0(full_dir, "/output_", df_type, ".dat"))
    # data contains columns Descendants (or AlleleCount in older versions) and pop_size
    colnames(temp)[colnames(temp) == "AlleleCount"] <- "Descendants"
    calc_VAF <- function(data){
      alpha <- 1
      coverage <- data$pop_size * alpha
      VAF <- data$Descendants / coverage
      return(VAF)
    }
    temp$pop_size <- (df_out %>% filter(Generation == max(Generation)) %>% select(NumCells))$NumCells
    temp$VAF <- calc_VAF(temp)
    if(!is.na(vaf_cut_off)) {
      temp <- temp %>% filter(VAF >= vaf_cut_off | Population > 0)
    }
    # add parameter columns:
    if(nrow(temp) == 0){
      temp <- NULL
    } else {
      temp <- cbind(df_pars, temp)
    }
  } else if (df_type %in% c("driver_phylo")){
    df_driver_phylo <- fread(file_driver_phylo)
    
    df_driver_phylo <- filter(df_driver_phylo, CellsPerSample == -1, NumSamples == 1, SampleDepth == -1)
    df_driver_phylo <- df_driver_phylo[!duplicated(df_driver_phylo), ]
    temp <- df_driver_phylo
    # add parameter columns:
    if(nrow(temp) == 0){
      temp <- NULL
    } else {
      temp <- cbind(df_pars, temp)
    }
  } else {
    stop("no valid df_type argument was passed")
  }
  
  # filter if requested:
  temp <- filter_by_generation_or_numcells(temp, full_dir, generation, numcells)
  
  print(paste0("Result of combine_dfs has dimensions ", dim(temp)[1], " x ", dim(temp)[2]), quote = FALSE)
  
  return(temp)
}

#' Create a composite dataframe for a batch of simulations, derived from multiple 
#' data files per simulation.
#' 
#' @param input_dir base input directory name
#' @param include_diversities boolean whether to include diversity metrics
#' @param df_type which dataframes to combine
#' @param vaf_cut_off exclude genotypes with vaf lower cut off from combined_df
#' @param generation Generation at which to filter (default NA corresponds to no filtering)
#' @param numcells Number of cells at which to filter (default NA corresponds to no filtering)
#' 
#' @return a combined dataframe
#' 
#' @importFrom data.table rbindlist
#' 
#' @export
#' 
#' @details If both \code{generation} and \code{numcells} are provided then \code{numcells} 
#' takes precedent.
#' 
#' @examples
#' all_output(system.file("example_batch", "", package = "demonanalysis", mustWork = TRUE))
#' all_output(system.file("example_batch", "", package = "demonanalysis", mustWork = TRUE), 
#' df_type = "driver_genotype_properties")
#' all_output(system.file("example_batch", "", package = "demonanalysis", mustWork = TRUE), 
#' df_type = "genotype_properties", vaf_cut_off = 0.002)
#' all_output(system.file("example_batch", "", package = "demonanalysis", mustWork = TRUE), 
#' df_type = "allele_counts", generation = 10)
all_output <- function(input_dir, include_diversities = TRUE, df_type = "output", vaf_cut_off = NA, generation = NA, numcells = NA) {
  
  df_type_list <- c("output", "driver_genotype_properties", "genotype_properties", 
                    "allele_counts", "driver_allele_counts", "genotype_counts", "driver_genotype_counts", "driver_phylo")
  stopifnot(df_type %in% df_type_list)
  
  pars_and_values <- parameter_names_and_values(input_dir)
  if(is.na(pars_and_values)[1]) stop("input_dir should contain results of a batch of simulations")
  pars <- pars_and_values$name
  final_values <- pars_and_values$final_value
  
  each_check <- function(x, res) {
    full_dir <- make_dir(input_dir, pars, x)
    msg <- final_error_message(full_dir)
    if(!identical(msg, character(0))) print(paste0(full_dir, " ", msg), quote = FALSE)
    else print("Failed")
  }
  apply_combinations(final_values, each_check)
  
  print("Finished checking", quote = FALSE)
  
  num_parameters <- count_parameters(input_dir)
  each_df <- function(x, res) {
    full_dir <- make_dir(input_dir, pars, x)
    msg <- final_error_message(full_dir)
    print(paste0(full_dir, " ", msg), quote = FALSE)
    if(!identical(msg, character(0))) if(msg == "Exit code 0") return(combine_dfs(full_dir, include_diversities, 
                                                                                  df_type, vaf_cut_off, generation, numcells, num_parameters))
    return(data.frame())
  }
  res <- rbindlist(apply_combinations(final_values, each_df))
  
  # report seed counts:
  print("Number of seeds:", quote = FALSE)
  print(count_seeds(res, num_parameters))
  
  return(res)
}

#' Count the number of rows in a dataframe per parameter set.
#' 
#' @param data dataframe
#' @param num_parameters number of parameters, accounting for the first set of columns in the dataframe
#' 
#' @return a dataframe listing number of rows per parameter set
#' 
#' @import dplyr
#' @export
#' 
#' @examples
#' count_seeds(sum_df, 16)
count_seeds <- function(data, num_parameters) {
  col_nums <- c(1:num_parameters)
  col_nums <- col_nums[col_nums != which(colnames(data) == "seed")]
  res <- data %>% group_by_at(col_nums) %>%
    summarise(num_seeds = length(unique(seed)))
  return(res$num_seeds)
}

#' Get summary metrics for one or more simulations
#' 
#' @param data dataframe
#' @param start_size_range vector of NumCells at time of initial measurement for forecasting
#' @param gap_range vector of projection periods (gap between time of initial measurement and second measurement)
#' @param final_size waiting time is measured until tumour reaches this NumCells value
#' @param num_parameters number of parameters, accounting for the first set of columns in the dataframe
#' 
#' @return a dataframe with one row for each unique combination of parameter values (including seed), 
#' "gap" and "start_size", (i.e. it summarises over time) and which has added columns "start_time" 
#' (proportional time until NumCells reached start_size), "end_time" (proportional time until NumCells 
#' reached end_size), "waiting_time" (difference between start_time and end_time), and "outcome" 
#' (NumCells after "gap")
#' 
#' @import dplyr
#' @export
#' 
#' @examples
#' num_parameters = count_parameters(system.file("extdata", "", 
#' package = "demonanalysis", mustWork = TRUE))
#' comb_df <- combine_dfs(system.file("extdata", "", package = "demonanalysis", mustWork = TRUE))
#' get_summary(comb_df, c(100, 300), c(0.35, 0.65), 1000, num_parameters)
get_summary <- function(data, start_size_range, gap_range, final_size, num_parameters) {
  summary <- data.frame()
  data <- data %>% group_by_at(1:num_parameters)
  for(start_size in start_size_range) {
    for(gap in gap_range) {
      if(start_size < final_size) {
        new_summary1 <- data %>% 
          filter(NumCells >= start_size, !is.na(DriverDiversity), !is.na(NumClones)) %>% 
          summarise(start_time = min(gen_adj, na.rm = TRUE))
        new_summary2 <- data %>% 
          filter(NumCells >= final_size, !is.na(DriverDiversity), !is.na(NumClones)) %>% 
          summarise(end_time = min(gen_adj, na.rm = TRUE))
        new_summary12 <- merge(new_summary1, new_summary2, all.x = TRUE)
        new_summary12 <- new_summary12 %>% 
          mutate(waiting_time = end_time - start_time)
      }
      else {
        new_summary12 <- data %>% 
          filter(NumCells >= start_size, !is.na(DriverDiversity), !is.na(NumClones)) %>% 
          summarise(waiting_time = NA, start_time = NA, end_time = NA)
      }
      new_summary3 <- data %>% 
        filter(NumCells > start_size, !is.na(DriverDiversity), !is.na(NumClones)) %>% 
        filter(gen_adj < min(gen_adj, na.rm = TRUE) + gap) %>% 
        summarise(outcome = max(NumCells, na.rm = TRUE))
      new_summary3a <- data %>% 
        filter(NumCells > start_size, !is.na(DriverDiversity), !is.na(NumClones)) %>% 
        summarise(outcome = max(NumCells, na.rm = TRUE))
      new_summary3$outcome <- ifelse(new_summary3$outcome == new_summary3a$outcome, NA, new_summary3$outcome)
      new_summary4 <- data %>% 
        filter(NumCells > start_size, !is.na(DriverDiversity), !is.na(NumClones)) %>% 
        filter(gen_adj == min(gen_adj, na.rm = TRUE)) %>%
        mutate(gap = gap, start_size = start_size)
      summary <- rbind(summary, merge(merge(new_summary12, new_summary3, all.x = TRUE), new_summary4, all.x = TRUE))
    }
  }
  summary <- summary %>% 
    mutate(Ratio = DriverEdgeDiversity / DriverDiversity)
  
  # check number of rows:
  count1 <- dim(summary)[1]
  count2 <- sum(count_seeds(summary, num_parameters)) * length(start_size_range) * length(gap_range)
  if(count1 != count2) stop(paste0("Row count (", count1, ") is not as expected (", count2, ")."))
  
  # report number of replicates per parameter set:
  print("Number of seeds:", quote = FALSE)
  print(count_seeds(summary, num_parameters))
  
  return(summary)
}

#' Generic function to find a correlation between two columns of a dataframe
#' 
#' @param summary dataframe
#' @param factor1 char name of first column in the summary dataframe
#' @param factor2 char name of second column in the summary dataframe
#' @param result_name char name to give the result
#' @param min_count minimum number of items in each column (otherwise result will be NA)
#' 
#' @return Correlation between the two columns (or NA if either factor1 or factor2 contains NA, 
#' or if all values of factor1 or factor2 are identical).
#' 
#' @import dplyr
#' @import lazyeval
#' @importFrom stats setNames
#' @importFrom stats var
#' @export
#' 
#' @examples
#' s1 <- data.frame(a = 1:3, b = 1:3 * (1 + rnorm(3) / 10))
#' find_correlations(s1, "a", "b", "c", 3)
find_correlations <- function(summary, factor1, factor2, result_name, min_count) {
  summary %>% 
    mutate_(variance = interp(~var(var1), var1 = as.name(factor2))) %>% 
    filter(variance > 0) %>% # to avoid warnings when all values of factor2 are identical
    mutate_(count1 = interp(~length(var1), var1 = as.name(factor1)), 
            count2 = interp(~length(var2), var2 = as.name(factor2))) %>% 
    filter(count1 >= min_count, count2 >= min_count) %>% 
    summarise_(temp_name = interp(~cor(var1, var2, method = "spearman"), var1 = as.name(factor1), var2 = as.name(factor2))) %>% 
    rename_(.dots = setNames("temp_name", paste0(result_name)))
}

#' Generate summary dataframe of correlations with "outcome"
#' 
#' @param summary dataframe including columns named "seed", "Generation", "start_size", "start_time" and "outcome"
#' @param col_names_list char vector of column names in the summary dataframe
#' @param num_parameters number of parameters, accounting for the first set of columns in the dataframe
#' @param min_count minimum number of items in each column (otherwise result will be NA)
#' 
#' @return Dataframe with one row for each unique combination of parameter values, gap and start_size 
#' (i.e. it summarises over "seed"), and including columns containing the correlations between "outcome" 
#' and each variable in col_names_list.
#' 
#' @import dplyr
#' @importFrom stats var
#' @export
#' 
#' @examples
#' get_cor_summary(sum_df, c("DriverDiversity", "DriverEdgeDiversity"), 16, min_count = 5)
get_cor_summary <- function(summary, col_names_list, num_parameters, min_count) {
  col_nums <- c(1:num_parameters, which(colnames(summary) == "gap"), which(colnames(summary) == "start_size"))
  col_nums <- col_nums[col_nums != which(colnames(summary) == "seed")]
  summary <- summary %>% 
    group_by_at(col_nums) %>% 
    filter(!is.na(outcome)) %>% 
    filter(var(outcome) > 0)
  cor_summary <- summary %>% 
    summarise(mean_start_time = mean(start_time), 
              mean_DriverDiversity = mean(DriverDiversity),
              mean_outcome = mean(outcome), 
              mean_autocor = mean(mean_autocor), 
              num_seeds = n())
  result_names_list <- paste0("Cor_", col_names_list)
  cor_summary_list <- list()
  for(i in 1:length(col_names_list)) cor_summary_list[[i]] <- find_correlations(summary, "outcome", col_names_list[i], result_names_list[i], min_count)
  for(i in 1:length(col_names_list)) cor_summary <- merge(cor_summary, cor_summary_list[[i]], all.x = TRUE)
  return(cor_summary)
}

#' Generate summary dataframe of correlations with "waiting_time"
#' 
#' @param summary dataframe including columns named "seed", "Generation", "start_time", "start_size", "gap" and "waiting_time"
#' @param col_names_list char vector of column names in the summary dataframe
#' @param num_parameters number of parameters, accounting for the first set of columns in the dataframe
#' @param min_count minimum number of items in each column (otherwise result will be NA)
#' 
#' @return Dataframe with one row for each unique combination of parameter values and start_size 
#' (i.e. it summarises over "seed"), and including columns containing the correlations between "waiting_time" 
#' and each variable in col_names_list.
#' 
#' @import dplyr
#' @importFrom stats var
#' @export
#' 
#' @examples
#' wait_cor_summary <- get_wait_cor_summary(sum_df, 
#' c("DriverDiversity", "DriverEdgeDiversity"), 16, min_count = 5)
#' depth_wait_cor_summary <- get_wait_cor_summary(sum_df, 
#' c(paste0("DriverDiversityFrom1SamplesAtDepth", 0:10), 
#' paste0("DriverDiversityFrom4SamplesAtDepth", 0:10)), 
#' 16, min_count = 5)
get_wait_cor_summary <- function(summary, col_names_list, num_parameters, min_count) {
  col_nums <- c(1:num_parameters, which(colnames(summary) == "start_size"))
  col_nums <- col_nums[col_nums != which(colnames(summary) == "seed")]
  summary <- summary %>% 
    group_by_at(col_nums) %>% 
    filter(gap == min(gap, na.rm = TRUE)) %>% # choice of gap value doesn't affect the result
    filter(!is.na(waiting_time)) %>% 
    filter(var(waiting_time) > 0)
  cor_summary <- summary %>% 
    summarise(mean_start_time = mean(start_time), 
              mean_DriverDiversity = mean(DriverDiversity), 
              mean_waiting_time = mean(waiting_time), 
              mean_autocor = mean(mean_autocor), 
              num_seeds = n())
  result_names_list <- paste0("Cor_", col_names_list)
  cor_summary_list <- list()
  for(i in 1:length(col_names_list)) cor_summary_list[[i]] <- find_correlations(summary, "waiting_time", col_names_list[i], result_names_list[i], min_count)
  for(i in 1:length(col_names_list)) cor_summary <- merge(cor_summary, cor_summary_list[[i]], all.x = TRUE)
  
  cor_summary <- arrange(cor_summary, K, migration_type, migration_edge_only, start_size)
  
  return(cor_summary)
}

#' Don't use this; it's just for refreshing data files on Rob's machine
#' 
#' @return NA
#' 
#' @importFrom readr read_delim
#' @export
refresh_data_files <- function() {
  driver_matrix <- read_delim("~/Documents/GitHub/demonanalysis/inst/extdata/driver_matrix.dat", "\t", col_names = FALSE)
  save(driver_matrix, file="data/driver_matrix.RData")
  driver_phylo <- read_delim("~/Documents/GitHub/demonanalysis/inst/extdata/driver_phylo.dat", "\t")
  save(driver_phylo, file="data/driver_phylo.RData")
  output_allele_cum_dist <- read_delim("~/Documents/GitHub/demonanalysis/inst/extdata/output_allele_cum_dist.dat", "\t")
  save(output_allele_cum_dist, file="data/output_allele_cum_dist.RData")
  output_allele_hist <- read_delim("~/Documents/GitHub/demonanalysis/inst/extdata/output_allele_hist.dat", "\t")
  save(output_allele_hist, file="data/output_allele_hist.RData")
  output_passengersgrid <- read_delim("~/Documents/GitHub/demonanalysis/inst/extdata/output_passengersgrid.dat", "\t", col_names = FALSE)
  save(output_passengersgrid, file="data/output_passengersgrid.RData")
  output <- read_delim("~/Documents/GitHub/demonanalysis/inst/extdata/output.dat", "\t")
  save(output, file="data/output.RData")
  parameters <- read_delim("~/Documents/GitHub/demonanalysis/inst/extdata/parameters.dat", "\t")
  save(parameters, file="data/parameters.RData")
  output_allele_counts <- read_delim(system.file("extdata", "output_allele_counts.dat", package = "demonanalysis", mustWork = TRUE), "\t")
  save(output_allele_counts, file="data/output_allele_counts.RData")
  output_driver_allele_counts <- read_delim(system.file("extdata", "output_driver_allele_counts.dat", package = "demonanalysis", mustWork = TRUE), "\t")
  save(output_driver_allele_counts, file="data/output_driver_allele_counts.RData")
  output_driver_genotype_counts <- read_delim(system.file("extdata", "output_driver_genotype_counts.dat", package = "demonanalysis", mustWork = TRUE), "\t")
  save(output_driver_genotype_counts, file="data/output_driver_genotype_counts.RData")
  output_genotype_counts <- read_delim(system.file("extdata", "output_genotype_counts.dat", package = "demonanalysis", mustWork = TRUE), "\t")
  save(output_genotype_counts, file="data/output_genotype_counts.RData")
}
