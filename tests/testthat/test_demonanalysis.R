context("All tests")

library(ggmuller)
library(dplyr)

test_that("sweep sequences are consistent for proportions and breaks", {
  phylo <- filter(driver_phylo, CellsPerSample == -1)
  pop_df <- get_population_df(phylo)
  lag_gens <- round(length(unique(pop_df$Generation))/6)
  
  expect_equal(sweep_sequence(pop_df, lag_type = "proportions", breaks = 6), 
               sweep_sequence(pop_df, lag_type = "generations", lag_gens = lag_gens))
})

test_that("apply_combinations works for equal and unequal sequences", {
  expect_equal(sort(apply_combinations(c(2, 3), mean)), sort(apply_combinations(c(3, 2), mean)))
  expect_length(apply_combinations(c(2, 3), mean), 3*4)
  expect_length(apply_combinations(c(2, 2), mean), 3*3)
})

test_that("ancestry result has correct dimensions and NA counts", {
  edges1 <- data.frame(Parent = rep(0, 2), Identity = 1:2)
  anc <- ancestry(edges1)
  expect_equal(dim(anc)[1], 3)
  expect_equal(dim(anc)[2], 3)
  expect_equal(sum(is.na(anc)), 1)
})

test_that("dominant results", {
  edges1 <- data.frame(Parent = rep(0, 2), Identity = 1:2)
  anc <- ancestry(edges1)
  pop_df <- data.frame(Generation = rep(0:2, each = 3), 
                       Identity = rep(0:2, times = 3), 
                       Population = c(1, 0, 0, 2, 2, 0, 3, 5, 1))
  expect_equal(dominant(anc, pop_df, 0.1), 2)
  expect_equal(dominant(anc, pop_df, 0.5), 1)
  expect_equal(dominant(anc, pop_df, 0.5, 0), 0)
  expect_equal(dominant(anc, pop_df, 0, 0), 2)
})

test_that("find_correlations perfect fit", {
  s1 <- data.frame(a = 1:3, b = 1:3)
  expect_equal(find_correlations(s1, "a", "b", "c", 3)$c, 1)
})






