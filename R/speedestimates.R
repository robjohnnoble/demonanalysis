#' Expected migration rate per deme
#' 
#' @param i number of cells in the deme
#' @param m migration rate per cell, relative to birth rate
#' @param r1 reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param migration_type 0 or 1 (see details)
#' @param d migration distance relative to 1/sqrt(K) (see details)
#' @param K deme carrying capacity (only required if d is not NA)
#' 
#' @return The migration rate.
#' 
#' @details If migration_type = 0 (meaning that migration is attempted only after a birth event) 
#' then the migration rate will be multiplied by a factor of r2 * r1. Otherwise no adjustment will be made.
#' 
#' If d and K are set then the calculation assumes that min(d / sqrt(K), 1) is the probability 
#' that a migration attempt will land outside the deme. 
#' 
#' Importantly, this function does not account for the chance that the migrating cell will 
#' land in a deme that already occupied.
#' 
#' @export
#' 
#' @examples
#' # if r1 * r2 = 1 then results are identical for migration_type = 0 and migration_type = 1: 
#' lambda(1, 0.1, 0.9, 1/0.9, 0)
#' lambda(1, 0.1, 0.9, 1/0.9, 1)
lambda <- function(i, m, r1, r2, migration_type = 0, d = NA, K = NA) {
  m <- min(m, 1) # m is a probability, so it cannot exceed 1
  res <- m * i
  if(!is.na(d)) {
    if(is.na(K)) stop("If d is set then K must also be set.")
    res <- res * min(d / sqrt(K), 1)
  }
  if(migration_type == 0) return(r1 * r2 * res)
  else if(migration_type == 1) return(res)
  else stop("Invalid migration_type")
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
#' @param m migration rate per cell, relative to birth rate
#' @param r1 birth rate of cells in the destination deme, relative to reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param migration_type 0 or 1 (see details)
#' @param death_rate death rate of migrating cells, relative to death rate of cells in destination deme
#' @param d migration distance relative to 1/sqrt(K); if NA then d is set to sqrt(K)
#' 
#' @return The rate of migration followed by survival.
#' 
#' @export
#' 
#' @details If migration_type = 0 (meaning that migration rate is correlated with birth rate) 
#' then the migration rate will be multiplied by a factor of r2 * r1. Otherwise no adjustment will be made.
#' 
#' @examples 
#' lambda_invasion(1, 2, 0.1, 1, 1.1, 0)
#' lambda_invasion(1, 2, 0.1, 1, 1.1, 1)
lambda_invasion <- function(i, K, m, r1, r2, migration_type = 0, death_rate = 1, d = NA) {
  res <- lambda(i, m, r1, r2, migration_type, d, K) * rho(1, K, r2)
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

#' Migration waiting time from fully occupied deme
#' 
#' Mean time until successful migration from a deme that is fully occupied by the focus cell type
#' 
#' @param K deme carrying capacity
#' @param m migration rate per cell, relative to birth rate
#' @param r1 birth rate of cells in the destination deme, relative to reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param migration_type 0 or 1 (see details)
#' @param d migration distance relative to 1/sqrt(K); if NA then d is set to sqrt(K)
#' 
#' @return The waiting time.
#' 
#' @export
#' 
#' @details Meaning of migration_type:
#' \itemize{
#'   \item 0: migration, such that a migration attempt occurs only following a birth event
#'   \item 1: migration, such that migration is independent of birth
#'   }
#' 
#' @examples 
#' time_migration(2, 0.1, 1, 1.1, 0)
#' time_migration(2, 0.1, 1, 1.1, 1)
time_migration <- function(K, m, r1, r2, migration_type, d = NA) {
  return(1 / lambda_invasion(K, K, m, r1, r2, migration_type, d = d))
}

#' Conditional fixation time
#' 
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

#' Moran growth waiting time
#' 
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

#' Migration waiting time (any population size)
#' 
#' Expected time of first successful migration, which may occur before fixation
#' 
#' @param r1 birth rate of cells in the destination deme, relative to reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param K deme carrying capacity
#' @param m migration rate per cell, relative to birth rate
#' @param migration_type 0 or 1 (see details)
#' @param d migration distance relative to 1/sqrt(K); if NA then d is set to sqrt(K)
#' 
#' @return The waiting time.
#' 
#' @export
#' 
#' @details Meaning of migration_type:
#' \itemize{
#'   \item 0: migration, such that a migration attempt occurs only following a birth event
#'   \item 1: migration, such that migration is independent of birth
#'   }
#' 
#' @examples 
#' time_expected(1, 1.1, 2, 0.1, 0)
#' time_expected(1, 1.1, 2, 0.1, 1)
#' 
#' # for comparison:
#' time_migration(2, 0.1, 1, 1.1, 0)
#' time_migration(2, 0.1, 1, 1.1, 1)
time_expected <- function(r1, r2, K, m, migration_type = 0, d = NA){
  
  if(K == 1) return(time_migration(K, m, r1, r2, migration_type, d))
  
  lambda_list <- sapply(1:K, lambda_invasion, K=K, m=m, r1=r1, r2=r2, migration_type = migration_type, d=d)
  t_list <- sapply(1:K, T_grow_j, r1=r1, r2=r2, K=K)
  
  l <- 2:K
  exponent_sum <- sapply(l, function(i) sum(lambda_list[1:(i-1)] * (t_list[1:(i-1)+1] - t_list[1:(i-1)])))
  term2 <- sum((1/lambda_list[l] - 1/lambda_list[l-1]) * exp(-exponent_sum[l-1]))

  return(1 / lambda_list[1] + term2)
}

#' Maximum limit on dispersal speed for migration model
#' 
#' @param K deme carrying capacity
#' @param m migration rate per cell, relative to birth rate
#' @param r1 birth rate of cells in the destination deme, relative to reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param migration_type 0 or 1 (see details)
#' @param symmetric whether migration occurs in both directions
#' @param two_dim whether to adjust for two-dimensional growth
#' @param d migration distance relative to 1/sqrt(K); if NA then d is set to sqrt(K)
#' 
#' @return The upper limit on the growth rate of the radius, measured in cells (not demes).
#' 
#' @export
#' 
#' @details Assumes that the time until fixation is negligible relative to the time until migration.
#' 
#' Meaning of migration_type:
#' \itemize{
#'   \item 0: migration, such that a migration attempt occurs only following a birth event
#'   \item 1: migration, such that migration is independent of birth
#'   }
#' 
#' @examples 
#' disp_rate_max(2, 0.1, 1, 1.1, 0)
#' disp_rate_max(2, 0.1, 1, 1.1, 1)
disp_rate_max <- function(K, m, r1, r2, migration_type = 0, symmetric = FALSE, two_dim = TRUE, d = NA) {
  m <- adjust_mig_rate(m, two_dim)
  if(!symmetric) return(sqrt(K) / time_migration(K, m, r1, r2, migration_type, d))
  else return((sqrt(K) / time_migration(K, m, r1, r2, migration_type, d) - sqrt(K) / time_migration(K, m, r1, 1 / r2, migration_type, d))/2)
}

#' Minimum limit on dispersal speed for migration model
#' 
#' @param r1 birth rate of cells in the destination deme, relative to reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param K deme carrying capacity
#' @param m migration rate per cell, relative to birth rate
#' @param migration_type 0 or 1 (see details)
#' @param symmetric whether migration occurs in both directions
#' @param two_dim whether to adjust for two-dimensional growth
#' @param d migration distance relative to 1/sqrt(K); if NA then d is set to sqrt(K)
#' 
#' @return The lower limit on the growth rate of the radius, measured in cells (not demes).
#' 
#' @export
#' 
#' @details Assumes that a cell type cannot attempt to migrate until it reaches fixation.
#' 
#' Meaning of migration_type:
#' \itemize{
#'   \item 0: migration, such that a migration attempt occurs only following a birth event
#'   \item 1: migration, such that migration is independent of birth
#'   }
#' 
#' @examples 
#' disp_rate_min(1, 1.1, 2, 0.1, 0)
#' disp_rate_min(1, 1.1, 2, 0.1, 1)
disp_rate_min <- function(r1, r2, K, m, migration_type = 0, symmetric = FALSE, two_dim = TRUE, d = NA) {
  m <- adjust_mig_rate(m, two_dim)
  if(!symmetric) return(sqrt(K) / (time_fixation(r1, r2, K) + time_migration(K, m, r1, r2, migration_type, d)))
  else return((sqrt(K) / (time_fixation(r1, r2, K) + time_migration(K, m, r1, r2, migration_type, d)) - sqrt(K) / 
                 (time_fixation(r1, 1 / r2, K) + time_migration(K, m, r1, 1 / r2, migration_type, d)))/2)
}

#' Expected disperal speed for migration model
#' 
#' Expected disperal speed for migration model, allowing dispersal to occur before fixation
#' 
#' @param r1 birth rate of cells in the destination deme, relative to reference birth rate
#' @param r2 birth rate of migrating cells, relative to r1
#' @param K deme carrying capacity
#' @param m migration rate per cell, relative to birth rate
#' @param migration_type 0 or 1 (see details)
#' @param symmetric whether migration occurs in both directions
#' @param two_dim whether to adjust for two-dimensional growth
#' @param d migration distance relative to 1/sqrt(K)
#' 
#' @return The growth rate of the radius, measured in cells (not demes).
#' 
#' @export
#' 
#' @examples 
#' disp_rate(1, 1.1, 2, 0.1, 0)
#' disp_rate(1, 1.1, 2, 0.1, 1)
#' 
#' # no difference between migration_type = 0 and migration_type = 1 when r1*r2 = 1:
#' sapply(0:1, disp_rate, r1 = 0.5, r2 = 2, K = 32, m = 0.1)
#' 
#' # otherwise expect a difference:
#' sapply(0:1, disp_rate, r1 = 0.5, r2 = 3, K = 32, m = 0.1)
disp_rate <- function(r1, r2, K, m, migration_type = 0, symmetric = FALSE, two_dim = TRUE, d = NA) {
  m <- adjust_mig_rate(m, two_dim)
  if(!symmetric) return(sqrt(K) / time_expected(r1, r2, K, m, migration_type, d))
  else return((sqrt(K) / time_expected(r1, r2, K, m, migration_type, d) - sqrt(K) / time_expected(r1, 1 / r2, K, m, migration_type, d))/2)
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
#' @return The growth rate of the radius, measured in cells (not demes).
#' 
#' @export
#' 
#' @details Assumes that a deme cannot attempt fission until its population size reaches K.
#' 
#' Meaning of migration_type:
#' \itemize{
#'   \item 2: fission, such that a fission attempt occurs only following a birth event
#'   \item 3: fission, such that fission is independent of birth
#'   }
#' 
#' @examples 
#' disp_rate_fission(1, 2, 0.1, NumCells = 1e6)
#' disp_rate_fission(1, 2, 0.1, migration_type = 3, NumCells = 1e6)
#' 
#' # compare observed and predicted dispersal rates:
#' df <- combine_dfs(system.file("extdata", "", package = "demonanalysis", mustWork = TRUE), 
#' include_diversities = FALSE)
#' rate_obs <- median(df$RadiusGrowthRate[which(df$NumCells > 400)])
#' rate_pred <- disp_rate_fission(1, 1, 1, migration_edge_only = 1)
#' rate_pred / rate_obs
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
    m_per_deme <- min(m * (K + 1), 1) # fission probability is per cell, so need to multiply by population size, which will be K+1 following a birth event
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

#' Demon migration rate
#' 
#' Quadratic formula used to set migration rate in demon.cpp
#' 
#' @param K deme carrying capacity
#' @param migration_type 0, 1, 2 or 3 (see details)
#' @param migration_edge_only whether migration occurs at the edge only
#' @param migration_rate_scales_with_K whether to divide by sqrt(K)
#' 
#' @return The dispersal speed.
#' 
#' @export
#' 
#' @details Meaning of migration_type:
#' \itemize{
#'   \item 0: migration, such that a migration attempt occurs only following a birth event
#'   \item 1: migration, such that migration is independent of birth
#'   \item 2: fission, such that a fission attempt occurs only following a birth event
#'   \item 3: fission, such that fission is independent of birth
#'   }
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

#' Two-dimensional adjusted migration rate
#' 
#' Adjust migration rate for probability that migrating cell lands in an already-occupied deme
#' 
#' @param m migration rate per cell, relative to birth rate
#' @param two_dim whether growth is in two dimensions
#' 
#' @return The migration rate.
#' 
#' @export
#' 
#' @details Adjustment for two-dimensional growth entails multiplying the migration rate by 
#' factors estimated from simulation results. First factor: the number of cells at the edge 
#' of the population, relative to the number expected if the population were a disc. 
#' Second factor: for cells at the edge, the average proportion of nearest neighbours that are empty.
#' 
#' @examples 
#' adjust_mig_rate(1, TRUE)
adjust_mig_rate <- function(m, two_dim) {
  if(two_dim) return(1.6 * 0.38 * m)
  else return(m)
}

#' Demon dispersal rate
#' 
#' A convenient wrapper for calculating dispersal rate using demon.cpp default parameter values
#' 
#' @param K deme carrying capacity
#' @param r2 birth rate of migrating cells
#' @param filled_grid whether the model is of a filled grid or an expanding tumour
#' @param migration_type 0, 1, 2, 3 (see details)
#' @param migration_edge_only whether migration occurs at the edge only
#' @param migration_rate_scales_with_K whether to divide migration rate by sqrt(K)
#' @param NumCells total population size (all demes); required only if migration_type = 2 or 3 and migration_edge_only = 0
#' 
#' @return The growth rate of the radius, measured in cells (not demes).
#' 
#' @export
#' 
#' @details Meaning of migration_type:
#' \itemize{
#'   \item 0: migration, such that a migration attempt occurs only following a birth event
#'   \item 1: migration, such that migration is independent of birth
#'   \item 2: fission, such that a fission attempt occurs only following a birth event
#'   \item 3: fission, such that fission is independent of birth
#'   }
#' 
#' If filled_grid = 1 then migration_type is set to 0.
#' 
#' @examples 
#' # compare migration_edge_only = 0 versus migration_edge_only = 1:
#' sapply(0:1, disp_rate_demon, K = 32, r2 = 1/0.9, migration_type = 0)
#' sapply(0:1, disp_rate_demon, K = 32, r2 = 1, migration_type = 2, NumCells = 1e6)
#' 
#' # filled grid example:
#' sapply(1:4, disp_rate_demon, r2 = 1.1, filled_grid = 1)
#' sapply(1:4, disp_rate_demon, r2 = 1.1, filled_grid = 1, migration_rate_scales_with_K = 0)
disp_rate_demon <- function(K, r2, filled_grid = 0, migration_type = 0, migration_edge_only = 0, migration_rate_scales_with_K = 1, NumCells = NA) {
  if(filled_grid) {
    m <- 1
    if(migration_rate_scales_with_K) m <- m / sqrt(K)
    return(disp_rate(1, r2, K, m, migration_type = 0, symmetric = FALSE, two_dim = TRUE, d = NA))
  }
  if(migration_type == 0 || migration_type == 1) {
    m <- mig_rate(K, migration_type, migration_edge_only, migration_rate_scales_with_K)
    return(disp_rate(0.9, r2, K, m, migration_type, symmetric = FALSE, two_dim = TRUE, d = NA))
  } else if(migration_type == 2 || migration_type == 3) {
    m <- mig_rate(K, migration_type, migration_edge_only, migration_rate_scales_with_K)
    return(disp_rate_fission(r2, K, m, migration_type, migration_edge_only, two_dim = TRUE, NumCells))
  } else {
    stop("Invalid type.")
  }
}
