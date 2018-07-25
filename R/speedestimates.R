#' Expected migration rate per deme
#' 
#' @param i number of cells in the deme
#' @param K deme carrying capacity
#' @param d migration distance relative to 1/sqrt(K)
#' @param m migration rate per cell, relative to birth rate
#' @param r1 reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param migration_type 0 or 1 (see details)
#' 
#' @return The migration rate.
#' 
#' @details Assumes that min(d / sqrt(K), 1) is the probability that a migration attempt will land outside the deme. 
#' Importantly, does not account for the chance that the migrating cell will land in an already-occupied deme.
#' 
#' If migration_type = 0 (meaning that migration rate is correlated with birth rate) 
#' then the migration rate will be multiplied by a factor of r2 * r1. Otherwise no adjustment will be made.
#' 
#' @export
#' 
#' @examples 
#' lambda(1, 2, 1, 0.1, 0.9, 1/0.9, 0)
#' lambda(1, 2, 1, 0.1, 0.9, 1/0.9, 1)
lambda <- function(i, K, d, m, r1, r2, migration_type = 0) {
  res <- m * i * min(d / sqrt(K), 1)
  if(migration_type == 0) return(r1 * r2 * res)
  else return(res)
}

#' Probability of fixation in a Moran process
#' 
#' @param i number of cells of the focus type
#' @param K deme carrying capacity
#' @param r relative fitness of the focus type
#' 
#' @return The fixation probability.
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
#' @param migration_type 0 or 1 (see details)
#' @param death_rate death rate of migrating cells, relative to death rate of cells in destination deme
#' 
#' @return The rate of migration followed by survival.
#' 
#' @export
#' 
#' @details If migration_type = 0 (meaning that migration rate is correlated with birth rate) 
#' then the migration rate will be multiplied by a factor of r2 * r1. Otherwise no adjustment will be made.
#' 
#' @examples 
#' lambda_invasion(1, 2, 1, 0.1, 1, 1.1, 0)
#' lambda_invasion(1, 2, 1, 0.1, 1, 1.1, 1)
lambda_invasion <- function(i, K, d, m, r1, r2, migration_type = 0, death_rate = 1) {
  res <- lambda(i, K, d, m, r1, r2, migration_type) * rho(1, K, r2)
  # note that the destination deme temporarily has population size K+1
  return(res * K / (K + death_rate))
}

#' Moran process transition rates
#' 
#' @param i initial population of the focus type
#' @param j new population of the focus type (should be i, i-1 or i+1)
#' @param r relative fitness of the focus type
#' @param K total population size
#' 
#' @return The transition rate.
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
#' @param migration_type 0, 1 or 4 (see details)
#' 
#' @return The waiting time.
#' 
#' @export
#' 
#' @details Meaning of migration_type:
#' 0: migration, such that migration rate is correlated with birth rate, and death rate is uniform
#' 1: migration, such that migration rate is independent of birth rate, and death rate is uniform
#' 
#' @examples 
#' time_migration(2, 1, 0.1, 1, 1.1, 0)
#' time_migration(2, 1, 0.1, 1, 1.1, 4)
time_migration <- function(K, d, m, r1, r2, migration_type) {
  return(1 / lambda_invasion(K, K, d, m, r1, r2, migration_type))
}

#' Mean time until one cell achieves fixation, given that it reaches fixation
#' 
#' @param r1 birth rate of non-focus cells, relative to reference birth rate
#' @param r2 birth rate of focus cell type
#' @param K deme carrying capacity
#' @param r_powers_shifted optional look-up vector
#' 
#' @return The waiting time.
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
#' @return The waiting time.
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
#' @param migration_type 0, 1 or 4 (see details)
#' 
#' @return The waiting time.
#' 
#' @export
#' 
#' @details Meaning of migration_type:
#' 0: migration, such that migration rate is correlated with birth rate, and death rate is uniform
#' 1: migration, such that migration rate is independent of birth rate, and death rate is uniform
#' 
#' @examples 
#' time_expected(1, 1.1, 2, 1, 0.1, 0)
#' time_expected(1, 1.1, 2, 1, 0.1, 4)
time_expected <- function(r1, r2, K, d, m, migration_type = 0){
  
  if(K == 1) return(time_migration(K, d, m, r1, r2, migration_type))
  
  lambda_list <- sapply(1:K, lambda_invasion, K=K, d=d, m=m, r1=r1, r2=r2, migration_type = migration_type)
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

#' Maximum limit on dispersal speed for migration model
#' 
#' @param K deme carrying capacity
#' @param d migration distance relative to 1/sqrt(K)
#' @param m migration rate per cell, relative to birth rate
#' @param r1 birth rate of cells in the destination deme, relative to reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param migration_type 0, 1 or 4 (see details)
#' @param symmetric whether migration occurs in both directions
#' @param two_dim whether to adjust for two-dimensional growth
#' 
#' @return The dispersal speed.
#' 
#' @export
#' 
#' @details Assumes that the time until fixation is negligible relative to the time until migration.
#' 
#' Meaning of migration_type:
#' 0: migration, such that migration rate is correlated with birth rate, and death rate is uniform
#' 1: migration, such that migration rate is independent of birth rate, and death rate is uniform
#' 
#' @examples 
#' disp_rate_max(2, 1, 0.1, 1, 1.1, 0)
#' disp_rate_max(2, 1, 0.1, 1, 1.1, 4)
disp_rate_max <- function(K, d, m, r1, r2, migration_type = 0, symmetric = FALSE, two_dim = TRUE) {
  m <- adjust_mig_rate(m, two_dim)
  if(!symmetric) return(sqrt(K) / time_migration(K, d, m, r1, r2, migration_type))
  else return((sqrt(K) / time_migration(K, d, m, r1, r2, migration_type) - sqrt(K) / time_migration(K, d, m, r1, 1 / r2, migration_type))/2)
}

#' Minimum limit on dispersal speed for migration model
#' 
#' @param r1 birth rate of cells in the destination deme, relative to reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param K deme carrying capacity
#' @param d migration distance relative to 1/sqrt(K)
#' @param m migration rate per cell, relative to birth rate
#' @param migration_type 0, 1 or 4 (see details)
#' @param symmetric whether migration occurs in both directions
#' @param two_dim whether to adjust for two-dimensional growth
#' 
#' @return The dispersal speed.
#' 
#' @export
#' 
#' @details Assumes that a cell type cannot attempt to migrate until it reaches fixation.
#' 
#' Meaning of migration_type:
#' 0: migration, such that migration rate is correlated with birth rate, and death rate is uniform
#' 1: migration, such that migration rate is independent of birth rate, and death rate is uniform
#' 
#' @examples 
#' disp_rate_min(1, 1.1, 2, 1, 0.1, 0)
#' disp_rate_min(1, 1.1, 2, 1, 0.1, 4)
disp_rate_min <- function(r1, r2, K, d, m, migration_type = 0, symmetric = FALSE, two_dim = TRUE) {
  m <- adjust_mig_rate(m, two_dim)
  if(!symmetric) return(sqrt(K) / (time_fixation(r1, r2, K) + time_migration(K, d, m, r1, r2, migration_type)))
  else return((sqrt(K) / (time_fixation(r1, r2, K) + time_migration(K, d, m, r1, r2, migration_type)) - sqrt(K) / 
                 (time_fixation(r1, 1 / r2, K) + time_migration(K, d, m, r1, 1 / r2, migration_type)))/2)
}

#' Expected disperal speed for migration model, allowing dispersal to occur before fixation
#' 
#' @param r1 birth rate of cells in the destination deme, relative to reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param K deme carrying capacity
#' @param d migration distance relative to 1/sqrt(K)
#' @param m migration rate per cell, relative to birth rate
#' @param migration_type 0, 1 or 4 (see details)
#' @param symmetric whether migration occurs in both directions
#' @param two_dim whether to adjust for two-dimensional growth
#' 
#' @return The dispersal speed.
#' 
#' @export
#' 
#' @examples 
#' disp_rate(1, 1.1, 2, 1, 0.1, 0)
#' disp_rate(1, 1.1, 2, 1, 0.1, 4)
disp_rate <- function(r1, r2, K, d, m, migration_type = 0, symmetric = FALSE, two_dim = TRUE) {
  m <- adjust_mig_rate(m, two_dim)
  if(!symmetric) return(sqrt(K) / time_expected(r1, r2, K, d, m, migration_type))
  else return((sqrt(K) / time_expected(r1, r2, K, d, m, migration_type) - sqrt(K) / time_expected(r1, 1 / r2, K, d, m, migration_type))/2)
}

#' Dispersal speed for deme fission model
#' 
#' @param r2 birth rate relative to reference birth rate
#' @param K deme carrying capacity
#' @param m migration rate per cell, relative to birth rate
#' @param migration_type 2 or 3 (see details)
#' @param migration_edge_only whether migration occurs at the edge only
#' @param two_dim whether to adjust for two-dimensional growth
#' @param NumCells total population size (all demes); required only if migration_type = 2 or 3 and migration_edge_only = 0
#' 
#' @return The dispersal speed.
#' 
#' @export
#' 
#' @details Assumes that a deme cannot attempt fission until its population size reaches K.
#' 
#' Meaning of migration_type:
#' 2: fission, such that fission rate is correlated with birth rate, and death rate is uniform
#' 3: fission, such that fission rate is independent of birth rate, and death rate is uniform
#' 
#' @examples 
#' disp_rate_fission(1, 2, 0.1, NumCells = 1e6)
#' disp_rate_fission(1, 2, 0.1, migration_type = 3, NumCells = 1e6)
#' 
#' \dontrun{
#' # find dispersal rate from simulation:
#' df <- combine_dfs(system.file("extdata", "", package = "demonanalysis", mustWork = TRUE))
#' num_parameters <- count_parameters(system.file("extdata", "", 
#' package = "demonanalysis", mustWork = TRUE))
#' df <- add_columns(df, num_parameters)
#' median(df$RadiusGrowthRate[which(df$NumCells > 400)])
#' 
#' # predicted dispersal rate is a bit too low in this case:
#' disp_rate_fission(1, 1, 1)
#' }
disp_rate_fission <- function(r2, K, m, migration_type = 2, migration_edge_only = 0, two_dim = TRUE, NumCells = NA) {
  # time_to_grow = time until deme population size reaches K
  if(K == 1) {
    time_to_grow <- 0 # time taken for population to grow from 1 to 1
  } else {
    time_to_grow <- 1 / r2 # time taken for population to grow from K / 2 to K
  }
  # fission_events_rate = number of attempted fission events per unit time:
  if(migration_type == 2) {
    birth_rate_per_deme <- r2 * K # in this case, fission events are proportional to birth events; deme population will usually be K
    m_per_deme <- min(m * (K + 1), 1) # fission rate is per cell, so need to multiply by population size, which will be K+1 following a birth event
    if(migration_edge_only) m_per_deme <- adjust_mig_rate(m_per_deme, two_dim)
    fission_events_rate <- birth_rate_per_deme * m_per_deme
  } else if(migration_type == 3) {
    if(migration_edge_only) m <- adjust_mig_rate(m, two_dim)
    if(K == 1) fission_events_rate <- 0 # fission can happen only when deme population > 1, which is rare in this case
    else fission_events_rate <- m * K # fission rate is per cell, so need to multiply by population size, which will usually be K
  } else {
    stop("Invalid migration_type")
  }
  time_to_migrate <- 1 / fission_events_rate
  if((migration_type == 2 || migration_type == 3) && migration_edge_only == 0) {
    if(is.na(NumCells)) stop("NumCells is missing, and is required in this case")
    return(1 / (time_to_grow + time_to_migrate) * 1/2 * sqrt(NumCells / pi))
  } else {
    return(sqrt(K) / (time_to_grow + time_to_migrate))
  }
}

#' Quadratic formula used to set migration rate in demon.cpp
#' 
#' @param K deme carrying capacity
#' @param migration_type 0, 1, 2, 3 or 4 (see details)
#' @param migration_edge_only whether migration occurs at the edge only
#' @param migration_rate_scales_with_K whether to divide by sqrt(K)
#' 
#' @return The dispersal speed.
#' 
#' @export
#' 
#' @details Meaning of migration_type:
#' 0: migration, such that migration rate is correlated with birth rate, and death rate is uniform
#' 1: migration, such that migration rate is independent of birth rate, and death rate is uniform
#' 2: fission, such that fission rate is correlated with birth rate, and death rate is uniform
#' 3: fission, such that fission rate is independent of birth rate, and death rate is uniform
#' 
#' @examples 
#' mig_rate(32)
#' mig_rate(32, migration_type = 2)
mig_rate <- function(K, migration_type = 0, migration_edge_only = 0, migration_rate_scales_with_K = 1) {
  if(migration_type == 0 || migration_type == 1) { # migration
    if(migration_edge_only == 0) {
      A <- -0.1593
      B <- -0.2868
      C <- 0.4646
    } else {
      A <- -0.2041
      B <- -0.14
      C <- 0.5761
    }
  } else if(migration_type == 2 || migration_type == 3) { # fission
    if(migration_edge_only == 0) {
      A <- -0.15
      B <- -0.8928
      C <- -2.3445
    } else {
      A <- -0.0431
      B <- -1.7951
      C <- 0.0726
    }
  } else {
    stop("Invalid migration_type")
  }
  m <- 10^(A * (log10(K))^2 + B * log10(K) + C)
  if(migration_rate_scales_with_K) m <- m / sqrt(K)
  if(migration_type == 3) m <- m * K
  return(m)
}

#' Adjust migration rate for probability that migrating cell lands in an already-occupied deme
#' 
#' @param m migration rate per cell, relative to birth rate
#' @param two_dim whether growth is in two dimensions
#' 
#' @return The migration rate.
#' 
#' @export
#' 
#' @details Adjustment for two-dimensional growth entails multiplying the migration rate by 0.4, which is
#' a factor estimated from simulation results.
#' 
#' @examples 
#' adjust_mig_rate(1, TRUE)
adjust_mig_rate <- function(m, two_dim) {
  if(two_dim) return(0.4 * m)
  else return(m)
}

#' A convenient wrapper for calculating dispersal rate for exanding tumours, using demon.cpp default parameter values
#' 
#' @param K deme carrying capacity
#' @param r2 birth rate of migrating cells
#' @param migration_type 0, 1, 2, 3 or 4 (see details)
#' @param migration_edge_only whether migration occurs at the edge only
#' @param NumCells total population size (all demes); required only if migration_type = 2 or 3 and migration_edge_only = 0
#' 
#' @return The dispersal speed.
#' 
#' @export
#' 
#' @details Assumes migration_type = 0, migration_rate_scales_with_K = 1.
#' 
#' @examples 
#' # compare migration_edge_only = 0 versus migration_edge_only = 1:
#' sapply(0:1, disp_rate_demon, K = 32, r2 = 1/0.9, migration_type = 0)
#' sapply(0:1, disp_rate_demon, K = 32, r2 = 1, migration_type = 2, NumCells = 1e6)
#' 
#' # no difference between migration_type = 0 and migration_type = 1 when r2 = 1/0.9 
#' # (because 0.9 is the birth rate of normal cells):
#' sapply(0:1, disp_rate_demon, K = 32, r2 = 1/0.9, migration_edge_only = 0)
#' sapply(0:1, disp_rate_demon, K = 32, r2 = 1/0.9, migration_edge_only = 1)
#' 
#' # otherwise expect a difference:
#' sapply(0:1, disp_rate_demon, K = 32, r2 = 2, migration_edge_only = 0)
#' sapply(0:1, disp_rate_demon, K = 32, r2 = 2, migration_edge_only = 1)
disp_rate_demon <- function(K, r2, migration_type = 0, migration_edge_only = 0, NumCells = NA) {
  if(migration_type == 0 || migration_type == 1) {
    return(disp_rate(0.9, r2, K, sqrt(K), mig_rate(K, migration_type, migration_edge_only), migration_type, symmetric = FALSE, two_dim = TRUE))
  } else if(migration_type == 2 || migration_type == 3) {
    return(disp_rate_fission(r2, K, mig_rate(K, migration_type, migration_edge_only), migration_type, migration_edge_only, two_dim = TRUE, NumCells))
  } else {
    stop("Invalid type.")
  }
}
