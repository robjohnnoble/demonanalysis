#' Apply a function to every combination of some sequences
#' 
#' @param vec vector of final values of the sequences (initial values are always zero)
#' @param fn function to apply to the values
#' @param ... other arguments passed to fn
#' 
#' @return result of applying fn to every combination of vec values
#' 
#' @export
#' 
#' @examples
#' apply_combinations(c(2, 3), mean)
apply_combinations <- function(vec, fn, ...){
  vecs <- mapply(seq, 0, vec, SIMPLIFY = FALSE)
  tmp <- do.call(expand.grid, vecs)
  apply(tmp, 1, fn, ...)
}

#' Attempt to read a tab-delimited file and return the contents, or NA if the file doesn't exist
#' 
#' @param file tab-delimited file
#' 
#' @return the contents, or NA if the file doesn't exist
#' 
#' @export
#' 
#' @examples
#' read_delim_special(system.file("extdata", "output_allele_counts.dat", 
#' package = "demonanalysis", mustWork = TRUE))
read_delim_special <- function(file) {
  if(file.exists(file)) out <- read_delim(file, "\t", trim_ws = TRUE)
  else out <- NA
  return(out)
}

#' Read a file containing grid states and process it into a dataframe for plotting.
#' 
#' @param file file name including path
#' @param trim how many rows and columns to remove; if trim < 0 (default) then all rows and columns containing NA are removed
#' 
#' @return a dataframe formatted for plotting
#' 
#' @export
#' @importFrom readr read_delim
#' 
#' @examples
#' image_df <- image_df_from_grid_file(system.file("extdata", 
#' "output_passengersgrid.dat", package = "demonanalysis", mustWork = TRUE))
image_df_from_grid_file <- function(file, trim = -1) {
  if(!file.exists(file)) {
    warning(paste0(file, " not found"))
    return(NA)
  }
  res <- read_delim(file, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  res <- data.matrix(res) # data frame to matrix
  res <- res[ , -ncol(res)] # drop extraneous columns
  # Mode <- function(x) {
  #   ux <- unique(x)
  #   ux[which.max(tabulate(match(x, ux)))]
  # }
  # background_val <- Mode(res[1, ]) # value denoting empty spaces
  if(length(res) > 1) {
    if(trim < 0) {
      res <- res[ , colSums(is.na(res)) < nrow(res)] # drop extraneous columns
      res <- res[rowSums(is.na(res)) < ncol(res), ] # drop extraneous rows
      # res <- res[ , colSums(res)/nrow(res) != background_val]
      # res <- res[rowSums(res)/ncol(res) != background_val , ]
    }
    else res <- res[(trim + 1):(nrow(res) - trim), (trim + 1):(ncol(res) - trim)]
  }
  df <- expand.grid(x = 1:ncol(as.data.frame(res)), y = 1:nrow(as.data.frame(res)))
  df$z <- as.vector(t(res))
  return(df)
}

#' Read a "phylo" dataframe and process it for plotting
#' 
#' @param file file name including path
#' 
#' @return a dataframe formatted for plotting
#' 
#' @export
#' @importFrom readr read_delim
#' @import dplyr
#' @import ggmuller
#' 
#' @examples
#' Muller_df <- muller_df_from_file(system.file("extdata", 
#' "driver_phylo.dat", package = "demonanalysis", mustWork = TRUE))
muller_df_from_file <- function(file) {
  if(!file.exists(file)) {
    warning(paste0(file, " not found"))
    return(NA)
  }
  phylo <- read_delim(file, "\t", escape_double = FALSE, trim_ws = TRUE)
  phylo <- filter(phylo, CellsPerSample == -1)
  edges <- get_edges(phylo)
  if(dim(edges)[1] == 0) edges <- NA
  pop_df <- get_population_df(phylo)
  pop_df <- pop_df %>% mutate(col_index = pop_df$Identity)
  pop_df$col_index[pop_df$col_index > 0] <- pop_df$col_index[pop_df$col_index > 0] %% 25 + 1
  pop_df$col_index <- as.character(pop_df$col_index)
  return(get_Muller_df(edges, pop_df))
}

#' Plot a grid from a properly formatted data frame.
#' 
#' @param image_df data frame formatted by image_df_from_grid_file
#' @param palette colour palette (default NA)
#' @param discrete whether to use a discrete or continuous colour scale (default FALSE)
#' @param add_legend whether to add a legend (default FALSE)
#' @param legend_title text for legend title (default none)
#' 
#' @return a plot object
#' 
#' @export
#' @import ggplot2
#' 
#' @examples
#' image_df <- image_df_from_grid_file(system.file("extdata", 
#' "output_passengersgrid.dat", package = "demonanalysis", mustWork = TRUE))
#' grid_plot(image_df)
grid_plot <- function(image_df, palette = NA, discrete = FALSE, add_legend = FALSE, legend_title = "") {
  if(length(image_df) == 1) return(ggplot(data.frame()) + theme_classic())
  h2 <- ggplot(image_df, aes(x, y, fill = z)) + 
    geom_raster() +
    theme(legend.position = ifelse(add_legend, "right", "none")) + 
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
  if(discrete) {
    if(!is.na(palette[1])) h2 <- h2 + scale_fill_manual(name = legend_title, values = palette) +
      scale_color_manual(values = palette)
  }
  else {
    if(is.na(palette[1])) {
      h2 <- h2 + scale_fill_distiller(name = legend_title, palette ="RdBu", direction = -1, na.value="white") + 
        scale_color_distiller(palette ="RdBu", na.value="white")
    }
    else {
      h2 <- h2 + scale_fill_distiller(name = legend_title, palette = palette, direction = -1) + 
        scale_color_distiller(palette = palette)
    }
  }
  return(h2)
}

#' Create a set of grid and Muller plots
#' 
#' @param path folder containing the input files
#' @param trim how many rows and columns to remove from grids; if trim < 0 (default) then all rows and columns containing NA are removed
#' @param output_dir folder in which to save the image file; if NA then plots are displayed on screen instead
#' @param output_filename name of output image file
#' @param file_type either "pdf" or "png" (other values default to "pdf")
#' 
#' @return either an image file or a plot displyed on screen
#' 
#' @export
#' @import ggmuller
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom grDevices png
#' 
#' @examples
#' plot_all_images(system.file("extdata", "", package = "demonanalysis", mustWork = TRUE))
plot_all_images <- function(path, output_filename = NA, file_type = "png", output_dir = NA, trim = -1) {
  if(substr(path, nchar(path), nchar(path)) != "/") path <- paste0(path, "/")
  if(!is.na(output_dir)) if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") output_dir <- paste0(output_dir, "/")
  
  Muller_df <- muller_df_from_file(paste0(path, "driver_phylo.dat"))
  if(class(Muller_df) != "data.frame") return(NA)
  
  long_palette <- c("#8A7C64", "#599861", "#89C5DA", "#DA5724", "#74D944", "#CE50CA", 
                    "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", 
                    "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                    "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", 
                    "#5E738F", "#D1A33D")
  dd <- (-1):25
  dd.col <- c("white", long_palette)
  names(dd.col)  <- dd
  
  h1 <- Muller_plot(Muller_df, colour_by = "col_index", palette = dd.col)
  h2 <- Muller_pop_plot(Muller_df, colour_by = "col_index", palette = dd.col)
  h3 <- Muller_plot(Muller_df, colour_by = "BirthRate", add_legend = TRUE)
  
  image_df <- image_df_from_grid_file(paste0(path, "output_driversgrid.dat"), trim)
  image_df[which(image_df$z > 0), "z"] <- as.character(image_df[which(image_df$z > 0), "z"] %% 25 + 1)
  g1 <- grid_plot(image_df, palette = dd.col, discrete = TRUE)
  
  image_df <- image_df_from_grid_file(paste0(path, "output_birthratesgrid.dat"), trim)
  g2 <- grid_plot(image_df, add_legend = TRUE, legend_title = "BirthRate")
  
  image_df <- image_df_from_grid_file(paste0(path, "output_passengersgrid.dat"), trim)
  g3 <- grid_plot(image_df, add_legend = TRUE, legend_title = "Passengers")
  
  image_df <- image_df_from_grid_file(paste0(path, "output_popgrid.dat"), trim)
  image_df[which(image_df$z == 0), "z"] <- NA
  g4 <- grid_plot(image_df, add_legend = TRUE, legend_title = "Population")
  
  if(!is.na(output_filename)) print(paste0("Created all plots for file ", output_filename), quote = FALSE)
  
  if(!is.na(output_filename) & !is.na(output_dir)) {
    if(file_type == "png") png(paste0(output_dir,output_filename,".png"), width = 1000, height = 1000, res = 100)
    else pdf(paste0(output_dir,output_filename,".pdf"), width = 10, height = 10)
  }
  lay <- rbind(c(1,1,2),
               c(3,3,3),
               c(4,4,5),
               c(NA,6,7))
  print(grid.arrange(h1, g1, h2, h3, g2, g3, g4, layout_matrix = lay, heights = c(1, 1, 0.75, 0.75)))
  if(!is.na(output_filename) & !is.na(output_dir)) dev.off()
  
  if(!is.na(output_filename)) print("Saved the plot", quote = FALSE)
}

#' Plot allele count versus origin time, coloured by birth rate
#' 
#' @param file file containing columns "AlleleCount", "OriginTime and "BirthRate"
#' @param log if TRUE then y-axis will be log-transformed (default FALSE)
#' 
#' @return plot displyed on screen
#' 
#' @export
#' 
#' @examples
#' plot_allelecount_vs_origintime(system.file("extdata", "output_genotype_properties.dat", 
#' package = "demonanalysis", mustWork = TRUE))
plot_allelecount_vs_origintime <- function(file, log = FALSE) {
  if(!file.exists(file)) {
    warning(paste0(file, " not found"))
    plot(0, type = 'n', axes = FALSE, ann = FALSE)
    return(NA)
  }
  df <- read_delim_special(file)
  
  q <- ggplot(filter(df, AlleleCount > 0), aes(OriginTime, AlleleCount, size = BirthRate, colour = log10(BirthRate))) + 
    geom_point(alpha = 0.5) + 
    scale_color_continuous(low = "blue", high = "red") + 
    theme_classic()
  
  if(log) q <- q + scale_y_continuous(trans = 'log10')
  
  print(q)
}

#' Plot counts of variant allele frequencies on linear scales
#' 
#' @param file file containing columns "Frequency" and "Count"
#' @param generation Generation at which to make the measurement (default NA corresponds to the final Generation)
#' @param ... other parameters passed to plot
#' 
#' @return plot displyed on screen
#' 
#' @export
#' 
#' @examples
#' plot_counts(system.file("extdata", "output_allele_counts.dat", 
#' package = "demonanalysis", mustWork = TRUE), ylim = c(0, 10))
plot_counts <- function(file, generation = NA, ...) {
  if(!file.exists(file)) {
    warning(paste0(file, " not found"))
    plot(0, type = 'n', axes = FALSE, ann = FALSE)
    return(NA)
  }
  df <- read_delim_special(file)
  if("Generation" %in% colnames(df)) {
    if(is.na(generation)) generation <- max(df$Generation)
    df <- filter(df, Generation == generation)
  }
  hist <- with(df, hist(rep(x = Frequency, times = Count), plot = FALSE, breaks = seq(0, 1, length = 100)))
  plot(hist, xlim = c(0, 1), ylab = "count", main = "", ...)
  abline(v = 0.1, lty = 2, col = "red")
}

#' Plot a histogram of variant allele frequencies with logit x-axis and log y-axis
#' 
#' @param file file containing columns "Frequency" and "Count"
#' @param generation Generation at which to make the measurement (default NA corresponds to the final Generation)
#' @param ... other parameters passed to plot
#' 
#' @return plot displyed on screen
#' 
#' @export
#' @import dplyr
#' @importFrom stats plogis
#' @importFrom stats qlogis
#' @importFrom graphics axis
#' @importFrom graphics lines
#' 
#' @examples
#' plot_logit_freq_dist(system.file("extdata", "output_allele_counts.dat", 
#' package = "demonanalysis", mustWork = TRUE))
plot_logit_freq_dist <- function(file, generation = NA, ...) {
  if(!file.exists(file)) {
    warning(paste0(file, " not found"))
    plot(0, type = 'n', axes = FALSE, ann = FALSE)
    return(NA)
  }
  df <- read_delim_special(file)
  if("Generation" %in% colnames(df)) {
    if(is.na(generation)) generation <- max(df$Generation)
    df <- filter(df, Generation == generation)
  }
  
  df <- filter(df, Frequency < 1, Frequency > plogis(-14))
  
  logit_breaks <- plogis(-14 + 0:100 * 26 / 100)
  hist <- with(df, hist(rep(x = Frequency, times = Count), plot = FALSE, breaks = logit_breaks))
  
  plot(log10(hist$density) ~ qlogis(hist$mids), 
       xaxt = "n", yaxt = "n", 
       xlim = c(qlogis(1E-6), qlogis(0.9999)), 
       ylim = c(-6, 6),
       ylab = "density", ...)
  
  xshort <- c(1E-6, 1E-4, 1E-2, 0.5, 0.99, 0.9999)
  axis(1, at = qlogis(xshort), labels = xshort)
  yshort <- -4:10
  axis(2, at = yshort, labels = 10^yshort)
  xlong <- plogis(seq(-8, -1, length = 100))
  lines(-1.5+log10(1/xlong^2) ~ qlogis(xlong), col = "black", lty = 2, lwd = 2)
  xlong <- plogis(seq(1, 8, length = 100))
  lines(-1.5+log10(1/((xlong-1)*log10(1-xlong))) ~ qlogis(xlong), col = "black", lty = 2, lwd = 2)
}

#' Plot cumulative density of variant allele frequencies versus inverse allele frequency
#' 
#' @param file file containing columns "Frequency" and "Count"
#' @param generation Generation at which to make the measurement (default NA corresponds to the final Generation)
#' @param ... other parameters passed to plot
#' 
#' @return plot displyed on screen
#' 
#' @export
#' 
#' @importFrom TeachingDemos subplot
#' 
#' @examples
#' plot_cum_dist(system.file("extdata", "output_allele_counts.dat", 
#' package = "demonanalysis", mustWork = TRUE))
plot_cum_dist <- function(file, generation = NA, ...) {
  if(!file.exists(file)) {
    warning(paste0(file, " not found"))
    plot(0, type = 'n', axes = FALSE, ann = FALSE)
    return(NA)
  }
  df1 <- read_delim_special(file)
  if("Generation" %in% colnames(df1)) {
    if(is.na(generation)) generation <- max(df1$Generation)
    df1 <- filter(df1, Generation == generation)
  }
  
  cum_count <- function(df, d) {
    indices <- which(df$Frequency > d)
    if(length(indices) == 0) return(0)
    else return(sum(df[indices, "Count"]))
  }
  
  InverseFrequency <- seq(1/0.5, 1/0.1, by = 0.1)
  CumulativeCount <- sapply(1/InverseFrequency, cum_count, df = df1)
  
  if(max(CumulativeCount) > 0) {
    min_y <- 0
    max_y <- max(CumulativeCount)
  }
  else {
    min_y <- 0
    max_y <- 1
  }
  
  plot(CumulativeCount ~ InverseFrequency, xaxt = "n", ylim = c(min_y, max_y), ylab = "cumulative count", ...)
  InverseFrequency_labels <- c(0.5, 0.25, 0.15, 0.1)
  axis(1, at = 1/InverseFrequency_labels, labels = paste0("1/", InverseFrequency_labels))
  
  # InverseFrequency <- seq(1, 10, by = 0.1)
  # CumulativeCount <- sapply(1/InverseFrequency, cum_count, df = df1)
  # 
  # subplot(plot(CumulativeCount ~ InverseFrequency, xlab = "", ylab = "", type = "l"), 
  #         x = 'topleft', inset = 0.05)
}

#' Get the first incomplete moment from frequency data
#' 
#' @param sizes for example, midpoints of a histogram
#' @param counts counts corresponding to the sizes
#' @param threshold lower bound of sizes to include
#' 
#' @return the first incomplete moment
#' 
#' @export
#' 
#' @examples
#' df_test <- data.frame(size = 1:20, count = exp(-(1:20)))
#' first_inc_moment(df_test$size, df_test$count, 2)
first_inc_moment <- function(sizes, counts, threshold) {
  mean_size <- sum(sizes * counts)
  sum1 <- sum(sizes[which(sizes >= threshold)] * counts[which(sizes >= threshold)])
  return(1 / mean_size * sum1)
}

#' Plot first incomplete moment from frequency data
#' 
#' @param sizes for example, midpoints of a histogram
#' @param counts counts corresponding to the sizes
#' @param max_size maximum size (default 1)
#' @param ... other parameters passed to plot
#' 
#' @return plot displyed on screen
#' 
#' @export
#' 
#' @examples
#' df_test <- data.frame(size = 1:20, count = exp(-(1:20)))
#' plot_first_inc_moment(df_test$size, df_test$count)
plot_first_inc_moment <- function(sizes, counts, max_size = 1, ...) { 
  if(length(sizes) <= 1) {
    plot(0, type = 'n', axes = FALSE, ann = FALSE)
    return(NA)
  }
  mom <- sapply(sizes, first_inc_moment, sizes = sizes, counts = counts)
  sizes <- sizes[which(mom > 0)]
  mom <- mom[which(mom > 0)]
  
  mom <- mom[which(sizes <= max_size)]
  sizes <- sizes[which(sizes <= max_size)]
  
  plot(mom ~ sizes, log = "y", ylab = "first incomplete moment", 
       xlim = c(0, max_size), ylim = c(min(mom), 1), ...)
}

#' Plot a set of charts representing allele frequencies and genotype sizes
#' 
#' @param path folder containing the input files
#' @param output_filename name of output image file
#' @param file_type either "pdf" or "png" (other values default to "pdf")
#' @param output_dir folder in which to save the image file; if NA then plots are displayed on screen instead
#' @param max_size maximum size (default NA corresponds to plotting frequencies, not sizes)
#' @param generation Generation at which to make the measurement (default NA corresponds to the final Generation)
#' 
#' @return plot displyed on screen
#' 
#' @export
#' 
#' @import dplyr
#' @importFrom readr read_delim
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom grDevices png
#' @importFrom graphics par
#' @importFrom graphics text
#' 
#' @examples
#' plot_all_charts(system.file("extdata", "", package = "demonanalysis", mustWork = TRUE))
plot_all_charts <- function(path, output_filename = NA, file_type = "png", output_dir = NA, max_size = NA, generation = NA) {
  if(substr(path, nchar(path), nchar(path)) != "/") path <- paste0(path, "/")
  if(!is.na(output_dir)) if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") output_dir <- paste0(output_dir, "/")
  
  if(!is.na(output_filename) & !is.na(output_dir)) {
    if(file_type == "png") png(paste0(output_dir,output_filename,".png"), width = 1100, height = 1100, res = 100)
    else pdf(paste0(output_dir,output_filename,".pdf"), width = 11, height = 11)
  }
  
  par(mfrow = c(4, 4))
  par(mgp = c(2.2, 1, 0))
  par(mar = c(3.8, 3.8, 0.8, 0.8))
  
  files_list <- c("output_allele_counts.dat", "output_driver_allele_counts.dat", "output_genotype_counts.dat", "output_driver_genotype_counts.dat")
  axis_lab <- c("mutation", "driver mutation", "genotype", "driver genotype")
  axis_lab2 <- c("clone", "driver clone", "genotype", "driver genotype")
  
  for(i in 1:4) {
    df1 <- read_delim_special(paste0(path, files_list[i]))
    if(is.na(generation)) generation <- max(df1$Generation)
    df1 <- filter(df1, Generation == generation)
    
    # plot 1:
    plot_counts(paste0(path, files_list[i]), xlab = paste0(axis_lab[i], " frequency"), ylim = c(0, 10))
    if(length(df1) > 1) div_alleles <- round(quadratic_diversity(df1[, c("Frequency", "Count")], 0.025, threshold = 0.1), 2)
    else div_alleles <- ""
    if(length(df1) > 1) text(1, 9, paste0("modes = ", div_alleles), pos = 2)
    
    # plot 2:
    plot_logit_freq_dist(paste0(path, files_list[i]), generation = generation, xlab = paste0(axis_lab[i], " frequency"))
    
    # plot 3:
    if(is.na(max_size)) plot_first_inc_moment(df1$Frequency, df1$Count, xlab = paste0(axis_lab[i], " frequency"))
    else plot_first_inc_moment(df1$Size, df1$Count, xlab = paste0(axis_lab2[i], " size"), max_size = max_size)
    
    # plot 4:
    plot_cum_dist(paste0(path, files_list[i]), generation = generation, xlab = paste0("inverse ", axis_lab[i], " frequency"))
  }
  
  if(!is.na(output_filename) & !is.na(output_dir)) dev.off()
}

#' Plot numers of cells with n mutations for each value of n
#' 
#' @param df data frame with columns "Generation", "CellsWith0Drivers" and "CellsWithXDrivers" for some values of X >= 1
#' 
#' @return plot displyed on screen
#' 
#' @export
#' 
#' @examples
#' plot_mutation_waves(output)
plot_mutation_waves <- function(df) {
  if(length(df) == 1) {
    plot(0, type = 'n', axes = FALSE, ann = FALSE)
    return(NA)
  }
  start_ind <- which(colnames(df) == "CellsWith1Drivers")
  end_ind <- dim(df)[2]
  plot(CellsWith0Drivers ~ Generation, data = df, type = "l", 
       ylim = c(10, max(df$NumCells)), log = "y", xlab = "Cell generations", ylab = "Number of cells")
  for(i in start_ind:end_ind) lines(df[ , i][[1]] ~ df$Generation)
}

#' Create a directory name including parameter names and values
#' 
#' @param input_dir base directory
#' @param pars vector of parameter names
#' @param indices vector of parameter values, of same length as pars
#' 
#' @return directory name
#' 
#' @export
#' 
#' @examples
#' make_dir("my_dir", c("migration_edge_only", "mu_driver_birth", "seed"), rep(0,3))
make_dir <- function(input_dir, pars, indices) {
  if(substr(input_dir, nchar(input_dir), nchar(input_dir)) == "/") input_dir <- substr(input_dir, 1, nchar(input_dir) - 1)
  
  if(length(indices) != length(pars)) stop("Unequal lengths of indices and pars.")
  for(i in 1:length(pars)) input_dir <- paste0(input_dir, "/", pars[i], "_", indices[i])
  input_dir <- paste0(input_dir, "/")
  return(input_dir)
}

#' Create an image file name including parameter names and values
#' 
#' @param prefix prefix for file name
#' @param pars vector of parameter names
#' @param indices vector of parameter values, of same length as pars
#' 
#' @return file name
#' 
#' @export
#' 
#' @examples
#' make_image_file_name("plot", c("migration_edge_only", "mu_driver_birth", "seed"), rep(0,3))
make_image_file_name <- function(prefix, pars, indices) {
  if(length(indices) != length(pars)) stop("Unequal lengths of indices and pars.")
  name <- prefix
  for(i in 1:length(pars)) name <- paste0(name, "_", pars[i], indices[i])
  return(name)
}

#' Return the final generation number of a simulation
#' 
#' @param input_dir directory name
#' 
#' @return generation number
#' 
#' @export
#' 
#' @examples
#' final_generation(system.file("extdata", "", package = "demonanalysis", mustWork = TRUE))
final_generation <- function(input_dir) {
  if(substr(input_dir, nchar(input_dir), nchar(input_dir)) == "/") input_dir <- substr(input_dir, 1, nchar(input_dir) - 1)
  
  res <- read_lines(paste0(input_dir, "/output.dat"))
  val <- strsplit(res[length(res)], "\t") # split the last line into a list of strings
  return(val[[1]][1])
}

#' Return the final line (or an earlier line) of an error log of a simulation
#' 
#' @param input_dir directory name
#' @param adjust number of lines prior to the final line (default 0)
#' 
#' @return line from error log file
#' 
#' @importFrom readr read_lines
#' @export
#' 
#' @examples
#' final_error_message(system.file("extdata", "", package = "demonanalysis", mustWork = TRUE))
#' final_error_message(system.file("extdata", "", package = "demonanalysis", mustWork = TRUE), 1)
final_error_message <- function(input_dir, adjust = 0) {
  if(substr(input_dir, nchar(input_dir), nchar(input_dir)) == "/") input_dir <- substr(input_dir, 1, nchar(input_dir) - 1)
  
  fname <- paste0(input_dir, "/error_log.dat")
  if(file.exists(fname))  {
    res <- read_lines(fname)
    return(res[length(res) - abs(adjust)])
  }
  else return("no file")
}

#' Return the final line (or an earlier line) of every error log in a batch of simulations.
#' 
#' @param input_dir base input directory name
#' @param adjust number of lines prior to the final line (default 0)
#' @param summary if TRUE then return a table instead of a list (default FALSE)
#' 
#' @return one line from each error log
#' 
#' @export
all_statuses <- function(input_dir, adjust = 0, summary = FALSE) {
  pars_and_values <- parameter_names_and_values(input_dir)
  if(is.na(pars_and_values)[1]) stop("input_dir should contain results of a batch of simulations")
  pars <- pars_and_values$name
  final_values <- pars_and_values$final_value
  
  each_msg <- function(x) {
    full_dir <- make_dir(input_dir, pars, x)
    return(final_error_message(full_dir, adjust))
  }
  stats <- apply_combinations(final_values, each_msg)
  
  if(!summary) return(stats)
  else {
    stats <- lapply(stats, function(x) if(identical(x, character(0))) "character(0)" else x)
    return(table(unlist(stats)))
  }
}

#' Create image files for every simulation in a batch
#' 
#' @param input_dir base input directory name
#' @param type what type of images to create: "plot", "chart" or "origintimes" (or a vector containing two or more of these strings)
#' @param file_type either "pdf" or "png" (other values default to "pdf")
#' @param output_dir folder in which to save the image files
#' @param max_size maximum size (default NA corresponds to plotting frequencies, not sizes)
#' @param generation Generation at which to make the measurement (default NA corresponds to the final Generation)
#' 
#' @return a set of image files
#' 
#' @export
create_plots_batch <- function(input_dir, type = "plot", file_type = "png", output_dir = NA, max_size = NA, generation = NA) {
  pars_and_values <- parameter_names_and_values(input_dir)
  if(is.na(pars_and_values)[1]) stop("input_dir should contain results of a batch of simulations")
  pars <- pars_and_values$name
  final_values <- pars_and_values$final_value
  
  each_plot <- function(x) {
    full_dir <- make_dir(input_dir, pars, x)
    msg <- final_error_message(full_dir)
    if(!identical(msg, character(0))) if(msg == "Exit code 0") {
      if("plot" %in% type) plot_all_images(full_dir, make_image_file_name("plot", pars, x), file_type, output_dir)
      if("chart" %in% type) plot_all_charts(full_dir, make_image_file_name("chart", pars, x), file_type, output_dir, max_size, generation)
      if("origintimes" %in% type) {
        if(!is.na(output_dir)) {
          fname <- make_image_file_name("origintimes", pars, x)
          if(file_type == "png") png(paste0(output_dir, "/", fname, ".png"), width = 700, height = 500, res = 100)
          else pdf(paste0(fname, ".pdf"), width = 7, height = 5)
        }
        plot_allelecount_vs_origintime(paste0(full_dir, "output_genotype_properties.dat"))
        if(!is.na(output_dir)) dev.off()
      }
    }
  }
  apply_combinations(final_values, each_plot)
}


