# outcome/waiting time versus diversity, coloured by K (JohnCode1)
# may need to cut out the data from after when tumours start to reach the max size
plot_outcome_versus_diversity <- function(cor_summary, output_filename = NA, file_type = "png", output_dir = NA) {
  if(!is.na(output_filename)) if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") output_dir <- paste0(output_dir, "/")
  
  if(!is.na(output_dir)) {
    if(file_type == "png") png(paste0(output_dir,output_filename,".png"), width = 1000, height = 1000, res = 100)
    else pdf(paste0(output_dir,output_filename,".pdf"), width = 500, height = 600)
  }
  
  start_size_range <- unique(cor_summary$start_size)
  K_range <- unique(cor_summary$K)
  
  cols <- rainbow(8)
  cols[5] <- "blue"
  
  par(mfrow=c(2, 3))
  par(mar=c(4.5, 4, 1.5, 1))
  for(start_size_val in start_size_range) {
    title <- paste("Measure at ", start_size_val, " cells", sep = "")
    for(K_val in K_range) {
      if(K_val == K_range[1]) {
        plot(1, type = "n", xlim = c(0, 1), ylim = c(-0.1, 1),
             main = title, xlab = "projection period", ylab = "correlation coefficient: diversity vs tumour size")
        if(start_size_val == start_size_range[1]) legend("topright", as.character(K_range), title = "K", 
                                                     xpd = TRUE, horiz = FALSE, inset = c(0.1, 0.03), bty = "n", lty = 1, col = cols)
      }
      lines(Cor_DriverDiversity ~ gap, data = filter(cor_summary, start_size == start_size_val, K == K_val), 
            col = cols[log2(K) + 1], lwd = 2)
      abline(h = 0, untf = FALSE, lty = 3)
    }
  }
  
  if(!is.na(output_dir)) dev.off()
}

# trajectories, coloured by K (JohnCode1)
plot_outcome_versus_waiting_time <- function(cor_summary, output_filename = NA, file_type = "png", output_dir = NA) {
  if(!is.na(output_filename)) if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") output_dir <- paste0(output_dir, "/")
  
  if(!is.na(output_dir)) {
    if(file_type == "png") png(paste0(output_dir,output_filename,".png"), width = 1000, height = 1000, res = 100)
    else pdf(paste0(output_dir,output_filename,".pdf"), width = 500, height = 600)
    }
  
  start_size_range <- unique(cor_summary$start_size)
  K_range <- unique(cor_summary$K)
  
  cols <- rainbow(8)
  cols[5] <- "blue"
  
  par(mfrow=c(1, 1))
  par(mar=c(4, 4, 2, 2))
  for(K_val in K_range) {
    if(K_val == K_range[1]) {
      plot(1, type = "n", xlim = c(0, 6E3), ylim = c(-1, 1),
           main = "", xlab = "tumour size at measurement", ylab = "correlation coefficient: diversity vs waiting time")
      legend("topright", as.character(K_range), title = "K", 
             xpd = TRUE, horiz = FALSE, inset = c(0.1, 0.03), bty = "n", lty = 1, col = cols)
      }
      lines(Cor_DriverDiversity ~ start_size, data = filter(cor_summary, K == K_val), 
            col = cols[log2(K) + 1], lwd = 2)
      abline(h = 0, untf = FALSE, lty = 3)
      }
  
  if(!is.na(output_dir)) dev.off()
}

# trajectories (relative time) coloured by K:
qplot(
  gen_adj,
  NumCells,
  group = interaction(K, seed),
  data = data,
  colour = factor(K),
  geom = "line")

# time adjusted trajectories by driver edge diversity (JohnCode1)
qp2 <- qplot(
  new_time,
  NumCells,
  group = interaction(seed, K, mig_rate),
  data = data,
  colour = div0,
  geom = "line")
qp2 + scale_colour_gradientn(colours=c(brewer.pal(5,"RdYlBu"), "#0000FF"), name = "diversity at periphery") +
  facet_grid(.~K, scales= "free") + 
  #scale_x_continuous(limits = c(0, 520), name = "time (years)", labels = 0:5) +
  scale_y_continuous(name = "tumour size") +
  theme_classic()

# correlations between waiting time and diversity measured at different depths, coloured by size at measurement (JohnCode1)


# example correlations of outcome/waiting time with edge diversity (JohnCode1)


# sweep metric against time (JohnCode2)


# violin plot of mean/skew autocorrelation (JohnCode2)


# growth trajectory coloured by mean sweep index for that simulation (JohnCode3; PlotCodeForJohn)


# R-squared by mean theta (JohnCode3)


# Mean rates (histograms_rates)






