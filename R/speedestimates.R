#' Expected migration rate per deme
#' 
#' @param i number of cells in the deme
#' @param K deme carrying capacity
#' @param d migration distance relative to 1/sqrt(K)
#' @param m migration rate per cell, relative to birth rate
#' @param r1 reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param selection if selection = "birth & migration" then the rate will be multiplied by a factor of r2 * r1
#' 
#' @return A number, corresponding to the rate.
#' 
#' @details Assumes that min(d / sqrt(K), 1) is the probability that a migration attempt will land outside the deme. 
#' Importantly, does not account for the chance that the migrating cell will land in an already-occupied deme.
#' 
#' @export
#' 
#' @examples 
#' lambda(1, 2, 1, 0.1, 1, 1.1)
lambda <- function(i, K, d, m, r1, r2, selection = "birth") {
  res <- m * i * min(d / sqrt(K), 1)
  if(selection == "birth & migration") return(r1 * r2 * res)
  else return(res)
}

#' Probability of fixation in a Moran process
#' 
#' @param i number of cells of the focus type
#' @param K deme carrying capacity
#' @param r relative fitness of the focus type
#' 
#' @return A number, corresponding to the probability.
#' 
#' @export
#' 
#' @examples 
#' rho(1, 8, 1.1)
rho <- function(i, K, r) {
  if(r == 1) return(i/K)
  if(r == 0) return(0)
  else return((1-1/r^i) / (1-1/r^K))
}

#' Rate of migration followed by survival
#' 
#' @param i number of cells of the focus type
#' @param K deme carrying capacity
#' @param d migration distance relative to 1/sqrt(K)
#' @param m migration rate per cell, relative to birth rate
#' @param r1 birth rate of cells in the destination deme, relative to reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param selection one of "birth", "death", "birth & migration"
#' 
#' @return A number, corresponding to the rate.
#' 
#' @export
#' 
#' @examples 
#' lambda_invasion(1, 2, 1, 0.1, 1, 1.1, "birth")
#' lambda_invasion(1, 2, 1, 0.1, 1, 1.1, "death")
lambda_invasion <- function(i, K, d, m, r1, r2, selection = "birth") {
  res <- lambda(i, K, d, m, r1, r2, selection) * rho(1, K, r2)
  # note that the destination deme temporarily has population size K+1
  if(selection == "death") return(res * K / (K + 1 / r2))
  else if(selection == "birth & migration" || selection == "birth")  return(res * K / (K + 1))
  else stop("Invalid selection option")
}

#' Moran process transition rates
#' 
#' @param i initial population of the focus type
#' @param j new population of the focus type (should be i, i-1 or i+1)
#' @param r relative fitness of the focus type
#' @param K total population size
#' 
#' @return A number, corresponding to the rate.
#' 
#' @export
#' 
#' @examples 
#' trans_rate(1, 2, 1.1, 4)
trans_rate <- function(i, j, r, K) {
  base <- (K-i)/(r*i+K-i)*i/K
  res <- NA
  if(j == 0) res <- 1
  if(i == K) res <- 1
  if(j == i-1) res <- base
  if(i == j) res <- 1 - base - r*base
  if(j == i+1) res <- r*base
  return(res)
}

#' Mean time until successful migration from a deme that is fully occupied by the focus cell type
#' 
#' @param K deme carrying capacity
#' @param d migration distance relative to 1/sqrt(K)
#' @param m migration rate per cell, relative to birth rate
#' @param r1 birth rate of cells in the destination deme, relative to reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param selection one of "birth", "death", "birth & migration"
#' 
#' @return A number, corresponding to the time.
#' 
#' @export
#' 
#' @examples 
#' time_migration(2, 1, 0.1, 1, 1.1, "birth")
#' time_migration(2, 1, 0.1, 1, 1.1, "death")
time_migration <- function(K, d, m, r1, r2, selection = "birth") {
  return(1 / lambda_invasion(K, K, d, m, r1, r2, selection))
}

#' Mean time until one cell achieves fixation, given that it reaches fixation
#' 
#' @param r1 birth rate of non-focus cells, relative to reference birth rate
#' @param r2 birth rate of focus cell type
#' @param K deme carrying capacity
#' @param r_powers_shifted optional look-up vector
#' 
#' @return A number, corresponding to the time.
#' 
#' @export
#' 
#' @details See Traulsen 2009, equation 1.38.
#' 
#' @examples 
#' time_fixation(1, 1.1, 4)
time_fixation <- function(r1, r2, K, r_powers_shifted = NA) {
  if(K == 1) return(0)
  if(r2 == 1) return(K - 1)
  
  if(is.na(r_powers_shifted[1])) r_powers_shifted <- r2^(0:K)
  
  l <- 1:(K-1)
  sum <- sum((1 - 1 / r_powers_shifted[l + 1]) * (r2 * l + K - l) * (1 - 1 / r_powers_shifted[K - l + 1]) / (l * (K - l)))
  if(r2^K > 1e20) {
    sum <- sum * 1 / (r2 - 1)
  } else {
    sum <- sum * r2^K / ((r2^K - 1) * (r2 - 1))
  }
  return(sum/r1) # relative to the doubling time of a cell with birth rate 1
}

#' Mean time from population size 1 until size j <= K in a Moran process
#' 
#' @param j final population size of the focus type
#' @param r1 birth rate of non-focus cells, relative to reference birth rate
#' @param r2 birth rate of focus cell type
#' @param K deme carrying capacity
#' 
#' @return A number, corresponding to the time.
#' 
#' @export
#' 
#' @details See Traulsen 2009, equations 1.38 and 1.39.
#' 
#' @examples 
#' T_grow_j(3, 1, 1.1, 4)
T_grow_j <- function(j, r1, r2, K) {
  if(j > K) return(Inf)
  if(j <= 1 || K <= 1) return(0)
  if(j == K) return(time_fixation(r1, r2, K))
  
  r_powers_shifted <- r2^(0:K) # need to add one to the index!
  t1 <- time_fixation(r1, r2, K, r_powers_shifted)
  
  if(r2 == 1) sum1 <- K / j - 1
  else {
    sum1 <- (r2^(1-j) - r2^(1-K)) / (r2 - 1) # sum of 1/r2^k from k=j to k=K-1
    sum1 <- sum1 * (1 - 1/r2) / (1 - 1/r2^j) # ratio rho(1)/rho(j)
  }
  
  l <- 1:(K-1)
  if(r2 == 1) sum2 <- K / j * sum((K - pmax(j, l)) / (K - l))
  else {
    sum2 <- sum((1 - 1 / r_powers_shifted[l + 1]) * (r2 * l + K - l) * (1 / r_powers_shifted[pmax(j, l) - l + 1] - 1 / r_powers_shifted[K - l + 1]) / (l * (K - l)))
    if(r2^j > 1e20) {
      sum2 <- sum2 / (r2 - 1)
    } else {
      sum2 <- sum2 * r2^j / ((r2^j - 1) * (r2 - 1))
    }
  }
  
  return((t1 * (1 + sum1) - sum2 / r1)) # relative to the doubling time of a cell with birth rate 1
}

#' Expected time of first successful migration, which may occur before fixation
#' 
#' @param r1 birth rate of cells in the destination deme, relative to reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param K deme carrying capacity
#' @param d migration distance relative to 1/sqrt(K)
#' @param m migration rate per cell, relative to birth rate
#' @param selection one of "birth", "death", "birth & migration"
#' 
#' @return A number, corresponding to the time.
#' 
#' @export
#' 
#' @examples 
#' time_expected(1, 1.1, 2, 1, 0.1, "birth")
#' time_expected(1, 1.1, 2, 1, 0.1, "death")
time_expected <- function(r1, r2, K, d, m, selection = "birth"){
  
  if(K == 1) return(time_migration(K, d, m, r1, r2, selection))
  
  lambda_list <- sapply(1:K, lambda_invasion, K=K, d=d, m=m, r1=r1, r2=r2, selection = selection)
  t_list <- sapply(1:K, T_grow_j, r1=r1, r2=r2, K=K)
  
  # first summand
  l <- 1:K
  part1 <- (lambda_list[l] * t_list[l] + 1) * exp(-lambda_list[l] * t_list[l])
  part2 <- (lambda_list[l] * t_list[l + 1] + 1) * exp(-lambda_list[l] * t_list[l + 1])
  exponent_sum <- sapply(1:(K-1), function(i) sum(lambda_list[1:(i-1)] * (t_list[1:(i-1)+1] - t_list[1:(i-1)])))
  exponent_sum[1] <- 0
  l <- 1:(K-1)
  term1 <- sum(1/lambda_list[l] * exp(-exponent_sum[l] + lambda_list[l] * t_list[l]) * (part1[l] - part2[l]))
  
  # second summand
  l <- 1:(K-1)
  exponent_sum <- sum(lambda_list[l] * (t_list[l+1] - t_list[l]))
  term2 <- 1/lambda_list[K] * exp(-exponent_sum) * (lambda_list[K] * t_list[K] + 1)
  
  return(term1 + term2)
}

#' Maximum limit on dispersal speed
#' 
#' @param K deme carrying capacity
#' @param d migration distance relative to 1/sqrt(K)
#' @param m migration rate per cell, relative to birth rate
#' @param r1 birth rate of cells in the destination deme, relative to reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param selection one of "birth", "death", "birth & migration"
#' @param symmetric whether migration occurs in both directions
#' 
#' @return A number, corresponding to the speed.
#' 
#' @export
#' 
#' @details Assumes that the time until fixation is negligible relative to the time until migration.
#' 
#' @examples 
#' disp_rate_max(2, 1, 0.1, 1, 1.1, "birth")
#' disp_rate_max(2, 1, 0.1, 1, 1.1, "death")
disp_rate_max <- function(K, d, m, r1, r2, selection = "birth", symmetric = FALSE) {
  if(!symmetric) return(sqrt(K) / time_migration(K, d, m, r1, r2, selection))
  else return((sqrt(K) / time_migration(K, d, m, r1, r2, selection) - sqrt(K) / time_migration(K, d, m, r1, 1 / r2, selection))/2)
}

#' Minimum limit on dispersal speed
#' 
#' @param r1 birth rate of cells in the destination deme, relative to reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param K deme carrying capacity
#' @param d migration distance relative to 1/sqrt(K)
#' @param m migration rate per cell, relative to birth rate
#' @param selection one of "birth", "death", "birth & migration"
#' @param symmetric whether migration occurs in both directions
#' 
#' @return A number, corresponding to the speed.
#' 
#' @export
#' 
#' @details Assumes that a cell type cannot attempt to migrate until it reaches fixation.
#' 
#' @examples 
#' disp_rate_min(1, 1.1, 2, 1, 0.1, "birth")
#' disp_rate_min(1, 1.1, 2, 1, 0.1, "death")
disp_rate_min <- function(r1, r2, K, d, m, selection = "birth", symmetric = FALSE) {
  if(!symmetric) return(sqrt(K) / (time_fixation(r1, r2, K) + time_migration(K, d, m, r1, r2, selection)))
  else return((sqrt(K) / (time_fixation(r1, r2, K) + time_migration(K, d, m, r1, r2, selection)) - sqrt(K) / (time_fixation(r1, 1 / r2, K) + time_migration(K, d, m, r1, 1 / r2, selection)))/2)
}

#' Expected disperal speed, allowing dispersal to occur before fixation
#' 
#' @param r1 birth rate of cells in the destination deme, relative to reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param K deme carrying capacity
#' @param d migration distance relative to 1/sqrt(K)
#' @param m migration rate per cell, relative to birth rate
#' @param selection one of "birth", "death", "birth & migration"
#' @param symmetric whether migration occurs in both directions
#' 
#' @return A number, corresponding to the speed.
#' 
#' @export
#' 
#' @examples 
#' disp_rate(1, 1.1, 2, 1, 0.1, "birth")
#' disp_rate(1, 1.1, 2, 1, 0.1, "death")
disp_rate <- function(r1, r2, K, d, m, selection = "birth", symmetric = FALSE) {
  if(!symmetric) return(sqrt(K) / time_expected(r1, r2, K, d, m, selection))
  else return((sqrt(K) / time_expected(r1, r2, K, d, m, selection) - sqrt(K) / time_expected(r1, 1 / r2, K, d, m, selection))/2)
}

#' Quadratic formula used to set migration rate in demon.cpp
#' 
#' @param K deme carrying capacity
#' @param edge_only whether migration occurs at the edge only
#' @param two_dim whether to adjust for two-dimensional growth
#' 
#' @return A number, corresponding to the speed.
#' 
#' @export
#' 
#' @details Assumes migration_type = 0, migration_edge_only = 0, migration_rate_scales_with_K = 1. 
#' Adjustment for two-dimensional growth entails multiplying the migration rate by 0.4, which is
#' a factor estimated from simulation results.
#' 
#' @examples 
#' mig_rate(8)
mig_rate <- function(K, edge_only = 0, two_dim = FALSE) {
  if(edge_only == 0) {
    A <- -0.1593
    B <- -0.2868
    C <- 0.4646
  } else {
    A <- -0.2041
    B <- -0.14
    C <- 0.5761
  }
  if(two_dim) {
    mult = 0.4 # approximately the right adjustment factor, based on simulation results
  } else {
    mult = 1
  }
  return(mult * 10^(A * (log10(K))^2 + B * log10(K) + C) / sqrt(K))
}

#' A convenient wrapper for calculating dispersal rate for parameter values used in demon.cpp
#' 
#' @param K deme carrying capacity
#' 
#' @param edge_only whether migration occurs at the edge only
#' @param two_dim whether to adjust for two-dimensional growth
#' @param r2 birth rate of migrating cells
#' @return A number, corresponding to the speed.
#' 
#' @export
#' 
#' @details Assumes migration_type = 0, migration_edge_only = 0, migration_rate_scales_with_K = 1.
#' 
#' @examples 
#' disp_rate_demon(8, 1/0.9)
disp_rate_demon <- function(K, r2, edge_only = 0, two_dim = FALSE) {
  disp_rate(1, r2, K, sqrt(K), mig_rate(K, edge_only, two_dim), "birth")
}
