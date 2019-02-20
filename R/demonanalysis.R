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
  vecs <- mapply(seq, 0, vec, SIMPLIFY = FALSE) # a list of n sequences, where n = length(vec)
  tmp <- do.call(expand.grid, vecs) # a data frame where each row is a permuation of values from the n sequences
  apply(tmp, 1, fn, ...) # the result of applying fn to each row of tmp
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
#' @param as_matrix whether to return a matrix or a dataframe
#' 
#' @return Either a dataframe formatted for plotting (as_matrix = FALSE) or a matrix (as_matrix = TRUE)
#' 
#' @export
#' @importFrom readr read_delim
#' 
#' @examples
#' image_df <- image_df_from_grid_file(system.file("extdata", 
#' "output_passengersgrid.dat", package = "demonanalysis", mustWork = TRUE))
image_df_from_grid_file <- function(file, trim = -1, as_matrix = FALSE) {
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
  
  if(as_matrix) return(res)
  
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

#' For cells at the boundary of a population, find mean proportion of nearest neighbour sites that are empty
#' 
#' @param mat matrix in which empty sites are denoted NA
#' 
#' @return a vector containing total number of empty neighbours, total number of edge cells, 
#' and average proportion of empty neighbours per edge cell
#' 
#' @export
#' 
#' @examples
#' mat <- image_df_from_grid_file(system.file("extdata", 
#' "output_popgrid.dat", package = "demonanalysis", mustWork = TRUE), as_matrix = TRUE)
#' prob_successful_migration(mat)
#' 
#' # Note that the edge of the grid counts as empty space:
#' mat <- matrix(c(NA,NA,NA,NA,
#' 1,1,1,1,
#' NA,NA,NA,NA), nrow = 3, byrow = TRUE)
#' prob_successful_migration(mat)
prob_successful_migration <- function(mat) {
  nrows <- dim(mat)[1]
  ncols <- dim(mat)[2]
  n_cells <- 0
  n_edge_cells <- 0
  n_empty_neighbours <- 0
  for(i in 1:nrows) for(j in 1:ncols) {
    temp <- 0
    if(!is.na(mat[i, j])) {
      if(i == 1 || is.na(mat[i - 1, j])) temp = temp + 1
      if(i == nrows || is.na(mat[i + 1, j])) temp = temp + 1
      if(j == 1 || is.na(mat[i, j - 1])) temp = temp + 1
      if(j == ncols || is.na(mat[i, j + 1])) temp = temp + 1
      n_cells <- n_cells + 1
    }
    n_empty_neighbours <- n_empty_neighbours + temp
    if(temp > 0) n_edge_cells <- n_edge_cells + 1
  }
  return(c(n_cells = n_cells, 
           n_empty_neighbours = n_empty_neighbours, 
           n_edge_cells = n_edge_cells, 
           prob = (n_empty_neighbours / n_edge_cells) / 4))
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
    if(!is.na(palette[1])) h2 <- h2 + scale_fill_manual(name = legend_title, values = palette, na.value="white") +
      scale_color_manual(values = palette, na.value="white")
  }
  else {
    if(is.na(palette[1])) {
      h2 <- h2 + scale_fill_distiller(name = legend_title, palette ="RdBu", direction = -1, na.value="white") + 
        scale_color_distiller(palette ="RdBu", na.value="white")
    }
    else {
      h2 <- h2 + scale_fill_distiller(name = legend_title, palette = palette, direction = -1, na.value="white") + 
        scale_color_distiller(palette = palette, na.value="white")
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
  dd <- 0:25
  dd.col <- long_palette
  names(dd.col) <- dd
  
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
#' @param file file containing columns "Descendants" (or "AlleleCount" in older versions), "OriginTime and "BirthRate"
#' @param log if TRUE then y-axis will be log-transformed (default FALSE)
#' @param colour_by character containing name of column by which to colour the plot
#' @param palette either a brewer palette or a vector of colours (if colour_by is categorical)
#' @param discrete whether to use a discrete or continuous colour scale (default FALSE)
#' 
#' @return plot displyed on screen
#' 
#' @importFrom grDevices col2rgb
#' 
#' @export
#' 
#' @examples
#' plot_allelecount_vs_origintime(system.file("extdata", "output_genotype_properties.dat", 
#' package = "demonanalysis", mustWork = TRUE))
#' plot_allelecount_vs_origintime(system.file("extdata", "output_genotype_properties.dat", 
#' package = "demonanalysis", mustWork = TRUE), colour_by = "DriverMutations", 
#' palette = c("black", "blue", "grey", "red"), discrete = TRUE)
plot_allelecount_vs_origintime <- function(file, log = FALSE, colour_by = "BirthRate", palette = NA, discrete = FALSE) {
  if(!file.exists(file)) {
    warning(paste0(file, " not found"))
    plot(0, type = 'n', axes = FALSE, ann = FALSE)
    return(NA)
  }
  
  df <- read_delim_special(file)
  df <- as.data.frame(df)
  colnames(df)[colnames(df) == "AlleleCount"] <- "Descendants"
  
  if(discrete) df[ , colour_by] <- as.factor(df[ , colour_by])
  
  direction <- 1
  if(is.na(palette[1])) {
    if(!discrete) {
      palette <- "RdBu"
      direction <- -1
    }
    else {
      long_palette <- c("#8A7C64", "#599861", "#89C5DA", "#DA5724", "#74D944", "#CE50CA", 
                        "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", 
                        "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                        "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", 
                        "#5E738F", "#D1A33D")
      palette <- rep(long_palette, ceiling(length(unique(df[ , colour_by])) / length(long_palette)))
    }
  }
  # test whether palette is a vector of colours; if not then we'll assume it's the name of a predefined palette:
  palette_named <- !min(sapply(palette, function(X) tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)))
  
  q <- ggplot(filter(df, Descendants > 0), aes_string("OriginTime", "Descendants", size = "BirthRate", colour = colour_by)) + 
    geom_point(alpha = 0.5) + 
    theme_classic()
  
  if(!discrete) {
    q <- q + 
      scale_color_distiller(palette = palette, direction = direction)
  }
  else {
    if(palette_named) {
      q <- q + 
        scale_color_brewer(palette = palette)
    }
    else {
      id_list <- sort(unique(select(df, colour_by))[[1]]) # list of legend entries, omitting NA
      q <- q + 
        scale_color_manual(values = palette)
    }
  }
  
  if(log) q <- q + scale_y_continuous(trans = 'log10')
  
  print(q)
}

#' Alternative to base hist function (using dplyr)
#' 
#' @param x a vector of values for which the histogram is desired
#' @param breaks a vector giving the breakpoints between histogram bins
#' @param counts optional vector of counts for each x value
#' 
#' @return data frame with counts and densities
#' 
#' @export
#' 
#' @examples
#' freq <- seq(0, 1, length = 20)
#' num <- rbinom(20, 10, 0.5)
#' breaks <- seq(0, 1, length.out = 5)
#' hist2(freq, breaks, num)
#' hist2(freq, breaks)
#' 
#' # equivalent using standard function:
#' hist(rep(x = freq, times = num), breaks, plot = FALSE)
#' hist(freq, breaks, plot = FALSE)
hist2 <- function(x, breaks, counts = 1) {
  bin_nums <- 1:(length(breaks) - 1)
  widths <- breaks - lag(breaks, 1)
  widths <- widths[!is.na(widths)]
  mids <- breaks[bin_nums] + widths / 2
  
  df <- data.frame(x = x, counts = counts)
  
  hist <- df %>% mutate(bin = cut(x, breaks = breaks, labels = bin_nums, include.lowest = TRUE)) %>%
    group_by(bin) %>% summarise(counts = sum(as.numeric(counts))) %>% 
    mutate(mids = mids[bin], density = counts / (sum(as.numeric(counts)) * widths[bin]))
  return(hist)
}

#' Plot counts of variant allele frequencies on linear scales
#' 
#' @param file_or_dataframe file or data frame containing columns "Frequency" and "Count"
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
plot_counts <- function(file_or_dataframe, generation = NA, ...) {
  if("data.frame" %in% class(file_or_dataframe)) df <- file_or_dataframe
  else {
    if(!file.exists(file_or_dataframe)) {
      warning(paste0(file_or_dataframe, " not found"))
      plot(0, type = 'n', axes = FALSE, ann = FALSE)
      return(NA)
    }
    df <- read_delim_special(file_or_dataframe)
  }
  if("Generation" %in% colnames(df)) {
    if(is.na(generation)) generation <- max(df$Generation)
    if(generation != "nofilter") df <- filter_by_generation_or_numcells(df, NA, generation, NA)
  }
  
  breaks <- seq(0, 1, length = 101)
  hist <- hist2(df$Frequency, breaks, df$Count)
  plot(hist$counts ~ hist$mids, type = "h", xlim = c(0, 1), ylab = "count", main = "", ...)
}

#' Plot a histogram of variant allele frequencies with logit x-axis and log y-axis
#' 
#' @param file_or_dataframe file or data frame containing columns "Frequency" and "Count"
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
plot_logit_freq_dist <- function(file_or_dataframe, generation = NA, ...) {
  if("data.frame" %in% class(file_or_dataframe)) df <- file_or_dataframe
  else {
    if(!file.exists(file_or_dataframe)) {
      warning(paste0(file_or_dataframe, " not found"))
      plot(0, type = 'n', axes = FALSE, ann = FALSE)
      return(NA)
    }
    df <- read_delim_special(file_or_dataframe)
  }
  if("Generation" %in% colnames(df)) {
    if(is.na(generation)) generation <- max(df$Generation)
    if(generation != "nofilter") df <- filter_by_generation_or_numcells(df, NA, generation, NA)
  }
  
  df <- filter(df, Frequency < 1, Frequency > plogis(-14))
  
  logit_breaks <- plogis(-14 + 0:100 * 26 / 100)
  
  hist <- hist2(df$Frequency, logit_breaks, df$Count)
  
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
#' @param file_or_dataframe file or data frame containing columns "Frequency" and "Count"
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
plot_cum_dist <- function(file_or_dataframe, generation = NA, max_y = NA, ...) {
  if("data.frame" %in% class(file_or_dataframe)) df1 <- file_or_dataframe
  else {
    if(!file.exists(file_or_dataframe)) {
      warning(paste0(file_or_dataframe, " not found"))
      plot(0, type = 'n', axes = FALSE, ann = FALSE)
      return(NA)
    }
    df1 <- read_delim_special(file_or_dataframe)
  }
  if("Generation" %in% colnames(df1)) {
    if(is.na(generation)) generation <- max(df1$Generation)
    if(generation != "nofilter") df1 <- filter_by_generation_or_numcells(df1, NA, generation, NA)
  }
  
  cum_count <- function(df, d) {
    indices <- which(df$Frequency > d)
    if(length(indices) == 0) return(0)
    else return(sum(df[indices, "Count"]))
  }
  
  InverseFrequency <- seq(1/0.5, 1/0.1, by = 0.1)
  CumulativeCount <- sapply(1/InverseFrequency, cum_count, df = df1)
  
  min_y <- 0
  if(is.na(max_y)){
    if(max(CumulativeCount) > 0) {
      max_y <- max(CumulativeCount)
    }
    else {
      max_y <- 1
    }    
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
#' @param condense either "discrete" or "continuous"
#' @param ... other parameters passed to plot
#' 
#' @return plot displyed on screen
#' 
#' @export
#' 
#' @examples
#' df_test <- data.frame(size = 1:20, count = exp(-(1:20)))
#' plot_first_inc_moment(df_test$size, df_test$count, max_size = 20)
plot_first_inc_moment <- function(sizes, counts, max_size = 1, condense = NA, ...) { 
  if(length(sizes) <= 1) {
    plot(0, type = 'n', axes = FALSE, ann = FALSE)
    return(NA)
  }
  
  if(!is.na(condense)) {
    if(condense == "discrete") {
      df <- data.frame(sizes = sizes, counts = counts)
      df <- group_by(df, sizes) %>%
        summarise(counts = sum(counts)) %>%
        ungroup()
      sizes <- df$sizes
      counts <- df$counts
    }
    else if(condense == "continuous") {
      breaks <- seq(0, max(sizes), length = 1e4 + 1)
      hist <- hist2(sizes, breaks, counts)
      sizes <- hist$mids
      counts <- hist$counts
    }
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
#' @param path_or_dflist folder containing the input files, or a list of data frames
#' @param output_filename name of output image file
#' @param file_type either "pdf" or "png" (other values default to "pdf")
#' @param output_dir folder in which to save the image file; if NA then plots are displayed on screen instead
#' @param max_size maximum size (default NA corresponds to plotting frequencies, not sizes)
#' @param generation Generation at which to make the measurement (default NA corresponds to the final Generation)
#' @param numcells Number of cells at which to make the measurement (default NA corresponds to the final size)
#' @param max_count Max value of y-axis in counts plot
#' @param num_parameters Number of parameters, accounting for the first set of columns in the dataframe
#' 
#' @return plot displyed on screen
#' 
#' @export
#' 
#' @details If both \code{generation} and \code{numcells} are provided then \code{numcells} 
#' takes precedent. A value for \code{num_parameters} is required only if the input data represents 
#' multiple simulations and either \code{generation} or \code{numcells} is specified.
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
#' plot_all_charts(list(output_allele_counts, output_driver_allele_counts, 
#' output_genotype_counts, output_driver_genotype_counts))
#' 
#' # combining results from multiple simulations:
#' df1 <- all_output(system.file("example_batch", "", package = "demonanalysis", mustWork = TRUE), 
#' df_type = "allele_counts", generation = 10)
#' df2 <- all_output(system.file("example_batch", "", package = "demonanalysis", mustWork = TRUE), 
#' df_type = "driver_allele_counts", generation = 10)
#' df3 <- all_output(system.file("example_batch", "", package = "demonanalysis", mustWork = TRUE), 
#' df_type = "genotype_counts", generation = 10)
#' df4 <- all_output(system.file("example_batch", "", package = "demonanalysis", mustWork = TRUE), 
#' df_type = "driver_genotype_counts", generation = 10)
#' num_parameters <- count_parameters(system.file("example_batch", "", 
#' package = "demonanalysis", mustWork = TRUE))
#' plot_all_charts(list(df1, df2, df3, df4), num_parameters = num_parameters)
plot_all_charts <- function(path_or_dflist, output_filename = NA, file_type = "png", output_dir = NA, max_size = NA, 
                            generation = NA, numcells = NA, max_count = 10, num_parameters = NA) {
  if("list" %in% class(path_or_dflist)) {
    input_list <- path_or_dflist
    path <- NA
  }
  else {
    path <- path_or_dflist
    if(substr(path, nchar(path), nchar(path)) != "/") path <- paste0(path, "/")
    input_list <- list("output_allele_counts.dat", "output_driver_allele_counts.dat", "output_genotype_counts.dat", "output_driver_genotype_counts.dat")
    input_list <- paste0(path, input_list)
  }
  
  if(!is.na(output_dir)) if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") output_dir <- paste0(output_dir, "/")
  
  if(!is.na(output_filename) & !is.na(output_dir)) {
    if(file_type == "png") png(paste0(output_dir,output_filename,".png"), width = 1100, height = 1100, res = 100)
    else pdf(paste0(output_dir,output_filename,".pdf"), width = 11, height = 11)
  }
  
  par(mfrow = c(4, 4))
  par(mgp = c(2.2, 1, 0))
  par(mar = c(3.8, 3.8, 0.8, 0.8))
  
  axis_lab <- c("mutation", "driver mutation", "genotype", "driver genotype")
  axis_lab2 <- c("clone", "driver clone", "genotype", "driver genotype")
  
  for(i in 1:4) {
    if("list" %in% class(path_or_dflist)) df1 <- input_list[[i]]
    else df1 <- read_delim_special(input_list[[i]])
    
    df1 <- filter_by_generation_or_numcells(df1, path = path, generation = generation, numcells = numcells, num_parameters = num_parameters)
    
    # plot 1:
    plot_counts(df1, xlab = paste0(axis_lab[i], " frequency"), generation = "nofilter", ylim = c(0, max_count))
    # if(length(df1) > 1) div_alleles <- round(quadratic_diversity(df1$Frequency, df1$Count, 0.025, threshold = 0.1), 2)
    # else div_alleles <- ""
    # if(length(df1) > 1) text(1, 0.9 * max_count, paste0("modes = ", div_alleles), pos = 2)
    
    # plot 2:
    plot_logit_freq_dist(df1, generation = "nofilter", xlab = paste0(axis_lab[i], " frequency"))
    
    # plot 3:
    if(is.na(max_size)) plot_first_inc_moment(df1$Frequency, df1$Count, xlab = paste0(axis_lab[i], " frequency"), condense = "continuous")
    else plot_first_inc_moment(df1$Size, df1$Count, xlab = paste0(axis_lab2[i], " size"), max_size = max_size, condense = "discrete")

    # plot 4:
    plot_cum_dist(df1, generation = "nofilter", xlab = paste0("inverse ", axis_lab[i], " frequency"))
  }
  
  if(!is.na(output_filename) & !is.na(output_dir)) dev.off()
}

#' Plot numbers of cells with n mutations for each value of n
#' 
#' @param df data frame with columns "Generation", "CellsWith0Drivers" and "CellsWithXDrivers" for some values of X >= 1
#' 
#' @return plot displyed on screen
#' 
#' @export
#' 
#' @examples
#' plot_mutation_waves(output)
plot_mutation_waves <- function (df) 
{
  if (length(df) == 1) {
    plot(0, type = "n", axes = FALSE, ann = FALSE)
    return(NA)
  }
  start_ind <- which(colnames(df) == "CellsWith1Drivers")
  end_ind <- dim(df)[2]
  plot(CellsWith0Drivers ~ Generation, data = df, type = "l", 
    ylim = c(10, max(df$NumCells)), log = "y", xlab = "Cell generations", 
    ylab = "Number of cells", col="brown")
  title("Mutation waves")
  for (i in start_ind:end_ind) lines(df[, i] ~ df$Generation, col = i)
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
  else return("no error_log.dat file")
}

#' Return the final line (or an earlier line) of every error log in a batch of simulations.
#' 
#' @param input_dir base input directory name
#' @param adjust number of lines prior to the final line (default 0)
#' @param summary if TRUE then return a table instead of a list (default FALSE)
#' @param with_names if TRUE (and if summary is FALSE) then include directory name with each result
#' 
#' @return one line from each error log
#' 
#' @export
all_statuses <- function(input_dir, adjust = 0, summary = FALSE, with_names = FALSE) {
  pars_and_values <- parameter_names_and_values(input_dir)
  if(is.na(pars_and_values)[1]) stop("input_dir should contain results of a batch of simulations")
  pars <- pars_and_values$name
  final_values <- pars_and_values$final_value
  
  each_msg <- function(x) {
    full_dir <- make_dir(input_dir, pars, x)
    if(!summary & with_names) return(paste0(full_dir, " ", final_error_message(full_dir, adjust)))
    else return(final_error_message(full_dir, adjust))
  }
  stats <- apply_combinations(final_values, each_msg)
  
  if(!summary) return(stats)
  else {
    stats <- lapply(stats, function(x) if(identical(x, character(0))) "So far no status" else x)
    if(length(stats) == 0) return("So far no statuses to report")
    return(table(unlist(stats)))
  }
}

#' Create image files for every simulation in a batch
#' 
#' @param input_dir base input directory name
#' @param type what type of images to create: "plot", "chart" or "origintimes" (or a vector containing two or more of these strings)
#' @param file_type either "pdf" or "png" (other values default to "pdf")
#' @param output_dir folder in which to save the image files
#' @param ... additional arguments passed to the plotting function
#' 
#' @return a set of image files
#' 
#' @export
create_plots_batch <- function(input_dir, type = "plot", file_type = "png", output_dir = NA, ...) {
  pars_and_values <- parameter_names_and_values(input_dir)
  if(is.na(pars_and_values)[1]) stop("input_dir should contain results of a batch of simulations")
  pars <- pars_and_values$name
  final_values <- pars_and_values$final_value
  
  each_plot <- function(x) {
    full_dir <- make_dir(input_dir, pars, x)
    msg <- final_error_message(full_dir)
    if(!identical(msg, character(0))) if(msg == "Exit code 0") {
      if("plot" %in% type) plot_all_images(full_dir, make_image_file_name("plot", pars, x), file_type, output_dir)
      if("chart" %in% type) plot_all_charts(full_dir, make_image_file_name("chart", pars, x), file_type, output_dir, ...)
      if("origintimes" %in% type) {
        if(!is.na(output_dir)) {
          fname <- make_image_file_name("origintimes", pars, x)
          if(file_type == "png") png(paste0(output_dir, "/", fname, ".png"), width = 700, height = 500, res = 100)
          else pdf(paste0(fname, ".pdf"), width = 7, height = 5)
        }
        plot_allelecount_vs_origintime(paste0(full_dir, "output_genotype_properties.dat"), ...)
        if(!is.na(output_dir)) dev.off()
      }
    }
  }
  apply_combinations(final_values, each_plot)
}

#' Create a data frame of genotype diversity versus time for a batch of simulations
#' 
#' @param input_dir base input directory name
#' 
#' @return a data frame
#' 
#' @export
#' 
#' @examples
#' driver_geno_div_df_batch(system.file("example_batch", "", 
#' package = "demonanalysis", mustWork = TRUE))
driver_geno_div_df_batch <- function(input_dir) {
  
  inv_Simpson_index <- function(p) 1 / sum(p*p)
  
  df <- all_output(input_dir, df_type = "driver_phylo")
  pars_and_values <- parameter_names_and_values(input_dir)
  par_names <- c("Generation", levels(pars_and_values$name))
  par_names[par_names == "log2_K"] <- "K"
  par_names[par_names == "log2_deme_carrying_capacity"] <- "K"
  
  sum_df <- group_by_at(df, par_names) %>% 
    mutate(Diversity = inv_Simpson_index(Population / sum(Population))) %>% 
    slice(1) %>% 
    ungroup()
  
  return(sum_df)
}

#' Plot genotype diversity versus time
#' 
#' @param sum_df data frame created by driver_geno_div_df_batch
#' @param output_filename name of output image file
#' @param file_type either "pdf" or "png" (other values default to "pdf")
#' @param output_dir folder in which to save the image file; if NA then plots are displayed on screen instead
#' @param log_y whether to apply a log10 transformation to the y axis (default FALSE)
#' @param facet1 column name for facet rows
#' @param facet2 column name for facet columns
#' @param height relative image height
#' @param width relative image width
#' @param ymax maximum of y axis
#' @param selected_seed seed number for which results will be highlighted
#' 
#' @return either an image file or a plot displyed on screen
#' 
#' @export
#' 
#' @details Lines are coloured by combinations of parameter values.
#' 
#' @examples
#' sum_df <- driver_geno_div_df_batch(system.file("example_batch", "", 
#' package = "demonanalysis", mustWork = TRUE))
#' plot_driver_geno_div(sum_df)
plot_driver_geno_div <- function(sum_df, output_filename = NA, file_type = "png", output_dir = NA, log_y = FALSE, 
                                       facet1 = "migration_type", facet2 = "migration_edge_only", height = 1, width = 1, ymax = 100, selected_seed = NA) {
  if(!is.na(output_dir)) if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") output_dir <- paste0(output_dir, "/")
  
  g1 <- ggplot(sum_df, aes(x = Generation, y = Diversity, group = interaction(K, migration_type, migration_edge_only, seed, s_driver_birth), colour = factor(K))) + 
    geom_line(col = "grey") + 
    facet_grid(reformulate(facet1, facet2)) + 
    theme_classic()
  
  if(!is.na(selected_seed)) g1 <- g1 + geom_line(data = filter(sum_df, seed == selected_seed), aes(x = Generation, y = Diversity), color = "black")
  
  if(log_y) g1 <- g1 + scale_y_log10(limits = c(1, ymax))
  else g1 <- g1 + scale_y_continuous(limits = c(1, ymax))
  
  if(!is.na(output_filename) & !is.na(output_dir)) {
    if(file_type == "png") png(paste0(output_dir,output_filename,".png"), width = 1000 * width, height = 1000 * height, res = 100)
    else pdf(paste0(output_dir,output_filename,".pdf"), width = 6 * width, height = 6 * height)
  }
  print(g1)
  if(!is.na(output_filename) & !is.na(output_dir)) dev.off()
  
  if(!is.na(output_filename)) print("Saved the plot", quote = FALSE)
}
