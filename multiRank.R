# Multi embryo rank orders

# Take the ranking code, and extend to arbitrary number of embryos

# What we want to know:
# How often would we get the lowest n aneuploid embryos from a total of k embryos
# with random aneuploidies and 0 dispersal?

library(tessera)
library(tidyverse)
library(parallel)
library(patchwork)
library(svglite)
library(data.table)

source("parameters.R")
source("functions.R")

# Assume dispersal is worst case - 0 for all
make.embryos <- function(aneu.vector, seed) {

  # Ensure sorted
  aneu.vector <- sort(aneu.vector)

  # Make list of embryos
  embryo.list <- lapply(aneu.vector, function(a) {
    tessera::Embryo(
      n.cells = 200,
      n.chrs = 1,
      prop.aneuploid = a,
      dispersal = 0,
      rng.seed = seed
    )
  })

  lapply(embryo.list, function(e) tessera::takeAllBiopsies(e, biopsy.size = 5, chromosome = 1))
}



# Sample from biopsies x times, and determine if the biopsies separate the n
# lowest aneuploidy embryos from the remainder
# The biopsy list, and the lowest n embryos we want to detect
# Return the accuracy scores of samples in which we separated the lowest n correctly
test.biopsy.list <- function(biopsy.list, x) {

  # preallocate matrix
  result <- matrix(nrow = x, ncol = length(biopsy.list))

  for (i in 1:x) {

    # Choose a random biopsy index from each embryo
    sampled.indexes <- sample(length(biopsy.list[[1]]), length(biopsy.list), replace = T)

    # Get the number of aneuploid cells in the chosen biopsy
    sampled.biopsies <- sapply(1:length(biopsy.list), function(i) biopsy.list[[i]][sampled.indexes[i]])

    # split into first n embryos and remainder
    # low.pool = sampled.biopsies[1:n]
    # high.pool = sampled.biopsies[(n+1):length(biopsy.list)]

    # Calculate the rank of each embryo from the biopsy
    result[i, ] <- rank(sampled.biopsies)
  }





  # We know embryos are sorted from lowest aneuploidy to highest, so we just
  # need to check if the first n embryo values are lower than all the remainder
  #
  # #  accuracy score: what fraction of the low pool are less than the high pool values?
  # sum(low.pool < min(high.pool))/length(low.pool)


  # return rank matrix
  result
}


# # this shows how frequently a set of embryos with aneuploidy levels given would be correctly separated into
# # two groups of low and high aneuploidy
# hist(unlist(result['mean',]), xlab = "Mean correct", main="Distribution for embryos 0, 0.2, 0.4, 0.8, 1")

# Given a vector of aneuploidies, test the separation
# Returns the mean fraction of biopsies that accurately
# separated the lowest two aneuploid embryos from the rest
test.aneuploidy.vector <- function(aneu.vector) {

  # Each seed in turn
  do.call(rbind, lapply(1:N.REPLICATES, function(s) {

    # Take all biopsies from pool of embryos
    biopsy.list <- make.embryos(aneu.vector, s)
    test.biopsy.list(biopsy.list = biopsy.list, x = 100)
  }))
}

calc.mean.ranks <- function(aneu.vector) {
  aneu.results <- test.aneuploidy.vector(aneu.vector)

  mean.result <- vector(mode = "numeric", length = length(aneu.vector))
  sd.result <- vector(mode = "numeric", length = length(aneu.vector))
  for (i in 1:length(aneu.vector)) {
    mean.result[i] <- mean(aneu.results[, i])
    sd.result[i] <- sd(aneu.results[, i])
  }
  data.frame("aneuploidy" = aneu.vector, "mean" = mean.result, "sd" = sd.result)
}


result1 <- calc.mean.ranks(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))

result2 <- calc.mean.ranks(c(0.2, 0.3, 0.4, 0.7, 0.7, 0.8, 0.8, 0.8, 0.8))

rank.plot.1 <- ggplot(result1, aes(x = rownames(result1), y = mean)) +
  geom_point() +
  geom_text(aes(label = aneuploidy * 100), nudge_y = 2.5) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd)) +
  scale_x_discrete() +
  scale_y_continuous(breaks = seq(0, nrow(result1), 1)) +
  labs(x = "Embryo in pool", y = "Mean rank from biopsies") +
  theme_classic()
# save.single.width(rank.plot.1, "Figure xxxx - Rank_multi_1", 85)

rank.plot.2 <- ggplot(result2, aes(x = rownames(result2), y = mean)) +
  geom_point() +
  geom_text(aes(label = aneuploidy * 100), nudge_y = 2.5) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd)) +
  scale_x_discrete() +
  scale_y_continuous(breaks = seq(0, nrow(result2), 1)) +
  labs(x = "Embryo in pool", y = "Mean rank from biopsies") +
  theme_classic()

save.double.width(rank.plot.1 + rank.plot.2 + plot_annotation(tag_levels = c("A")), "Figure xxxx - Rank_multi_2", 85)
