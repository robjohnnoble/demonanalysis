#' Read a file containing grid states and process it into a dataframe for plotting.
#' 
#' @param file file name including path
#' @param trim how many rows and columns to remove; if trim < 0 (default) then all rows and columns containing NA are removed
#' 
#' @return a dataframe formatted for plotting
#' 
#' @export
#' @import readr
#' 
#' @examples
#' image_df <- image_df_from_grid_file("data/output_passengersgrid.dat")
image_df_from_grid_file <- function(file, trim = -1) {
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
#' @import readr
#' @import dplyr
#' @import ggmuller
#' 
#' @examples
#' Muller_df <- muller_df_from_file("data/driver_phylo.dat")
muller_df_from_file <- function(file) {
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
#' @param palette colour palette
#' @param discrete whether to use a discrete or continuous colour scale
#' @param add_legend whether to add a legend
#' @param legend_title text for legend title
#' 
#' @return a plot object
#' 
#' @export
#' @import ggplot2
#' 
#' @examples
#' image_df <- image_df_from_grid_file("data/output_passengersgrid.dat")
#' grid_plot(image_df)
grid_plot <- function(image_df, palette, discrete = FALSE, add_legend = FALSE, legend_title = "") {
  h2 <- ggplot(image_df, aes(x, y, fill = z)) + 
    geom_raster() +
    theme(legend.position = ifelse(add_legend, "right", "none")) + 
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
  if(discrete) {
    h2 <- h2 + scale_fill_manual(name = legend_title, values = palette) +
      scale_color_manual(values = palette)
  }
  else {
    h2 <- h2 + scale_fill_distiller(name = legend_title, palette ="RdBu", direction = -1, na.value="white") + 
      scale_color_distiller(palette ="RdBu", na.value="white")
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
#' @import gridExtra
#' 
#' @examples
#' plot_all_images("data")
plot_all_images <- function(path, output_filename = "plot", file_type = "png", output_dir = NA, trim = -1) {
  if(substr(path, nchar(path), nchar(path)) != "/") path <- paste0(path, "/")
  if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") output_dir <- paste0(output_dir, "/")
  
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
  g4 <- grid_plot(image_df, add_legend = TRUE, legend_title = "Population")
  
  print("Created all plots")
  
  if(!is.na(output_dir)) {
    if(file_type == "png") png(paste0(output_dir,output_filename,".png"), width = 1000, height = 1000, res = 100)
    else pdf(paste0(output_dir,output_filename,".pdf"), width = 500, height = 600)
  }
  lay <- rbind(c(1,1,2),
               c(3,3,3),
               c(4,4,5),
               c(NA,6,7))
  print(grid.arrange(h1, g1, h2, h3, g2, g3, g4, layout_matrix = lay, heights = c(1, 1, 0.75, 0.75)))
  if(!is.na(output_dir)) dev.off()
  
  print("Saved the plot")
}

#' Plot a histogram of variant allele frequencies
#' 
#' @param df data frame containing "frequency" and "density" columns
#' 
#' @return plot displyed on screen
#' 
#' @export
#' 
#' @examples
#' library(readr)
#' output_allele_hist <- read_delim("data/output_allele_hist.dat", "\t", trim_ws = TRUE)
#' plot_allele_hist(output_allele_hist)
plot_allele_hist <- function(df) {
  plot(log10(density) ~ qlogis(frequency), data = df, 
       xaxt = "n", yaxt = "n", 
       xlim = c(qlogis(1E-5), qlogis(0.9999)), 
       ylim = c(-6, 6),
       xlab = "allele frequency", ylab = "density")
  xshort <- c(1E-4, 1E-2, 0.5, 0.99, 0.9999)
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
#' @param df data frame containing "inverse_frequency" and "cumulative_count" columns
#' 
#' @return plot displyed on screen
#' 
#' @export
#' 
#' @examples
#' library(readr)
#' output_allele_cum_dist <- read_delim("data/output_allele_cum_dist.dat", "\t", trim_ws = TRUE)
#' plot_allele_cum_dist(output_allele_cum_dist)
plot_allele_cum_dist <- function(df) {
  plot(cumulative_count ~ inverse_frequency, data = df,
       xlim = c(0, 100), ylim = c(0, 200), 
       xlab = "inverse allele frequency", ylab = "cumulative count")
}

#' Get a histogram of genotype sizes
#' 
#' @param file name of file containing lists of genotype sizes
#' 
#' @return histogram object
#' 
#' @export
#' 
#' @examples
#' hist1 <- get_genotype_sizes_hist("data/genotypes.dat")
get_genotype_sizes_hist <- function(file) {
  lastline <- function(filename) {
    out <- system(sprintf("wc -l %s", filename), intern = TRUE)
    n <- as.integer(sub(sprintf("[ ]*([0-9]+)[ ]%s", filename), "\\1", out))
    print(n)
    scan(filename, what="", skip = n - 1, nlines = 1, sep = "\n", quiet = TRUE)
  }
  geno_list <- lastline(file)
  geno_list <- strsplit(geno_list, "\t")
  geno_list <- as.numeric(geno_list[[1]])
  hist1 <- hist(geno_list, plot = FALSE, breaks = 50)
  return(hist1)
}

#' Plot a histogram of genotype sizes
#' 
#' @param hist histogram of genotype sizes
#' 
#' @return plot displyed on screen
#' 
#' @export
#' 
#' @examples
#' hist1 <- get_genotype_sizes_hist("data/genotypes.dat")
#' plot_genotype_sizes_hist(hist1)
plot_genotype_sizes_hist <- function(hist) {
  plot(hist$density / sum(hist$density) ~ hist$mids, log = "y", 
       xlim = c(0, 1E4), ylim = c(1E-5, 1), 
       xlab = "genotype size", ylab = "frequency")
}

#' Plot first incomplete moment of genotype sizes
#' 
#' @param hist histogram of genotype sizes
#' 
#' @return plot displyed on screen
#' 
#' @export
#' 
#' @examples
#' hist1 <- get_genotype_sizes_hist("data/genotypes.dat")
#' plot_first_inc_moment(hist1)
plot_first_inc_moment <- function(hist) { 
  first_inc_moment <- function(sizes, counts, n) {
    mean_size <- sum(sizes * counts)
    sum1 <- sum(sizes[which(sizes >= n)] * counts[which(sizes >= n)])
    return(1 / mean_size * sum1)
  }
  mom <- sapply(hist$mids, first_inc_moment, counts = hist$density, sizes = hist$mids)
  plot(mom ~ hist$mids, log = "y", 
       xlim = c(0, 1E4), ylim = c(1E-3, 1), 
       xlab = "genotype size", ylab = "first incomplete moment")
}

#' Plot a set of charts representing allele frequencies and genotype sizes
#' 
#' @param path folder containing the input files
#' @param output_dir folder in which to save the image file; if NA then plots are displayed on screen instead
#' @param output_filename name of output image file
#' @param file_type either "pdf" or "png" (other values default to "pdf")
#' 
#' @return plot displyed on screen
#' 
#' @export
#' 
#' @import readr
#' 
#' @examples
#' plot_all_charts("data")
plot_all_charts <- function(path, output_filename = "chart", file_type = "png", output_dir = NA) {
  if(substr(path, nchar(path), nchar(path)) != "/") path <- paste0(path, "/")
  if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") output_dir <- paste0(output_dir, "/")
  
  output_allele_hist <- read_delim(paste0(path, "output_allele_hist.dat"), "\t", trim_ws = TRUE)
  output_allele_cum_dist <- read_delim(paste0(path, "output_allele_cum_dist.dat"), "\t", trim_ws = TRUE)
  hist1 <- get_genotype_sizes_hist(paste0(path, "genotypes.dat"))
  
  if(!is.na(output_dir)) {
    if(file_type == "png") png(paste0(output_dir,output_filename,".png"), width = 700, height = 800, res = 100)
    else pdf(paste0(output_dir,output_filename,".pdf"), width = 7, height = 8)
  }
  
  par(mfrow = c(2, 2))
  
  plot_allele_hist(output_allele_hist)
  plot_allele_cum_dist(output_allele_cum_dist)
  plot_genotype_sizes_hist(hist1)
  plot_first_inc_moment(hist1)
  
  if(!is.na(output_dir)) dev.off()
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
#' library(readr)
#' output <- read_delim("data/output.dat", "\t", trim_ws = TRUE)
#' plot_mutation_waves(output)
plot_mutation_waves <- function(df) {
  start_ind <- which(colnames(df) == "CellsWith1Drivers")
  end_ind <- dim(df)[2]
  plot(CellsWith0Drivers ~ Generation, data = df, type = "l", 
       ylim = c(10, max(df$Population)), log = "y", xlab = "Cell generations", ylab = "Number of cells")
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
#' final_generation("data")
final_generation <- function(input_dir) {
  if(substr(input_dir, nchar(input_dir), nchar(input_dir)) == "/") input_dir <- substr(input_dir, 1, nchar(input_dir) - 1)
  
  res <- read_lines(paste0(input_dir, "/output.dat"))
  val <- strsplit(res[length(res)], "\t") # split the last line into a list of strings
  return(val[[1]][1])
}

#' Return the final line of an error log of a simulation
#' 
#' @param input_dir directory name
#' 
#' @return final line of error log
#' 
#' @export
#' 
#' @examples
#' final_error_message("data")
final_error_message <- function(input_dir) {
  if(substr(input_dir, nchar(input_dir), nchar(input_dir)) == "/") input_dir <- substr(input_dir, 1, nchar(input_dir) - 1)
  
  res <- read_lines(paste0(input_dir, "/error_log.dat"))
  return(res[length(res)])
}

#' Create image files for every simulation in a batch
#' 
#' @param input_dir base input directory name
#' @param output_dir output directory name
#' @param pars vector of parameter names
#' @param final_values vector of largest parameter values, of same length as pars
#' @param type what type of images to create: "plot" or "chart" or c("plot", "chart")
#' 
#' @return a set of image files
#' 
#' @export
#' 
#' @examples
#' 
create_plots_batch <- function(input_dir, output_dir, pars, final_values, type = "plot") {
  N <- length(pars)
  if(N != length(final_values)) stop("Unequal lengths of pars and final_values.")
  
  if(N == 1) for(a in 0:final_values[1]) {
    full_dir <- make_dir(input_dir, pars, a)
    msg <- final_error_message(full_dir)
    if(!identical(msg, character(0))) if(msg == "Exit code 0") {
      if("plot" %in% type) plot_all_images(full_dir, make_image_file_name("plot", pars, a), "png", output_dir)
      if("chart" %in% type) plot_all_charts(full_dir, make_image_file_name("chart", pars, a), "png", output_dir)
    }
  }
  
  if(N == 2) for(a in 0:final_values[1]) for(b in 0:final_values[2]) {
    full_dir <- make_dir(input_dir, pars, c(a, b))
    msg <- final_error_message(full_dir)
    if(!identical(msg, character(0))) if(msg == "Exit code 0") {
      if("plot" %in% type) plot_all_images(full_dir, make_image_file_name("plot", pars, c(a, b)), "png", output_dir)
      if("chart" %in% type) plot_all_charts(full_dir, make_image_file_name("chart", pars, c(a, b)), "png", output_dir)
    }
  }
  
  if(N == 3) for(a in 0:final_values[1]) for(b in 0:final_values[2]) 
    for(c in 0:final_values[3]) {
      full_dir <- make_dir(input_dir, pars, c(a, b, c))
      msg <- final_error_message(full_dir)
      if(!identical(msg, character(0))) if(msg == "Exit code 0") {
        if("plot" %in% type) plot_all_images(full_dir, make_image_file_name("plot", pars, c(a, b, c)), "png", output_dir)
        if("chart" %in% type) plot_all_charts(full_dir, make_image_file_name("chart", pars, c(a, b, c)), "png", output_dir)
      }
    }
  
  if(N == 4) for(a in 0:final_values[1]) for(b in 0:final_values[2]) 
    for(c in 0:final_values[3]) for(d in 0:final_values[4]) {
      full_dir <- make_dir(input_dir, pars, c(a, b, c, d))
      msg <- final_error_message(full_dir)
      if(!identical(msg, character(0))) if(msg == "Exit code 0") {
        if("plot" %in% type) plot_all_images(full_dir, make_image_file_name("plot", pars, c(a, b, c, d)), "png", output_dir)
        if("chart" %in% type) plot_all_charts(full_dir, make_image_file_name("chart", pars, c(a, b, c, d)), "png", output_dir)
      }
    }
  
  if(N == 5) for(a in 0:final_values[1]) for(b in 0:final_values[2]) 
    for(c in 0:final_values[3]) for(d in 0:final_values[4]) for(e in 0:final_values[5]) {
      full_dir <- make_dir(input_dir, pars, c(a, b, c, d, e))
      msg <- final_error_message(full_dir)
      if(!identical(msg, character(0))) if(msg == "Exit code 0") {
        if("plot" %in% type) plot_all_images(full_dir, make_image_file_name("plot", pars, c(a, b, c, d, e)), "png", output_dir)
        if("chart" %in% type) plot_all_charts(full_dir, make_image_file_name("chart", pars, c(a, b, c, d, e)), "png", output_dir)
      }
    }
}


