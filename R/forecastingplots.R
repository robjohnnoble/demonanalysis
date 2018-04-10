# outcome/waiting time versus diversity, coloured by K (JohnCode1)
plot_outcome_versus_diversity <- function(cor_summary, output_filename = NA, file_type = "png", output_dir = NA) {
  if(!is.na(output_filename)) if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") output_dir <- paste0(output_dir, "/")
  
  if(!is.na(output_dir)) {
    if(file_type == "png") png(paste0(output_dir,output_filename,".png"), width = 1000, height = 1000, res = 100)
    else pdf(paste0(output_dir,output_filename,".pdf"), width = 500, height = 600)
  }
  
  par(mfrow=c(2,3))
  par(mar=c(4.5,4,1.5,1))
  for(start_size in start_size_range) {
    plot_data <- cor_summary[which(cor_summary$start_size == start_size & cor_summary$mig_rate == 0.1), ]
    title <- paste("Measure at ", start_size, " cells", sep = "")
    for(K in K_range) {
      if(K == K_range[1]) {
        plot(Cor_Diversity ~ gap, data = plot_data, type = "l", col = "white", ylim = c(-0.1, 1),
             main = title, xlab = "projection period", ylab = "correlation coefficient: diversity vs tumour size")
        if(start_size == start_size_range[1]) legend("topright", as.character(K_range), title = "K", xpd = TRUE, horiz = FALSE, inset = c(0.1, 0.03), bty = "n", lty=1, col = cols)
      }
      lines(Cor_Diversity ~ gap, data = plot_data[which(plot_data$K == K), ], col = cols[log2(K)+1], lwd=2)
      abline(h = 0, untf = FALSE, lty = 3)
    }
  } # may need to cut out the data from after when tumours start to reach the max size
  #yellow is K2, blue is k16, pink is K128
  
  if(!is.na(output_dir)) dev.off()
}

# trajectories, coloured by K (JohnCode1)


# time adjusted trajectories by driver edge diversity (JohnCode1)


# correlations between waiting time and diversity measured at different depths, coloured by size at measurement (JohnCode1)


# example correlations of outcome/waiting time with edge diversity (JohnCode1)


# sweep metric against time (JohnCode2)


# violin plot of mean/skew autocorrelation (JohnCode2)


# growth trajectory coloured by mean sweep index for that simulation (JohnCode3; PlotCodeForJohn)


# R-squared by mean theta (JohnCode3)


# Mean rates (histograms_rates)






