# Test calculating rank order for mosaic biopsies

# given two embryos with 30% and 60% mosaicism (i.e. the midpoints of the "low"
# vs "high" mosaicism classes), what is the probability that they will be ranked
# in the correct order by the biopsies, and how does this depend on dispersal?

library(tessera)
library(tidyverse)
library(parallel)
library(patchwork)
library(svglite)
library(data.table)

source("parameters.R")
source("functions.R")

#' Calculate the ranks for a given combination of aneuploidy and dispersal
#'
#' @param p1 the fractional aneuploidy of embryo 1
#' @param p2 the fractional aneuploidy of embryo 2
#' @param d1 the dispersal of embryo 1
#' @param d2 the dispersal of embryo 2
#' @param seed1 the RNG seed for embryo 1
#' @param seed2 the RNG seed for embryo 2
#'
#' @return a named vector with the total biopsies in each rank class from all
#'  pairwise combinations of biopsies from the two embryos
calculate.ranks <- function(p1, p2, d1, d2, seed1, seed2) {
  low.embryo <- tessera::Embryo(
    n.cells = 200,
    n.chrs = 1,
    prop.aneuploid = p1,
    dispersal = d1,
    concordance = 1,
    rng.seed = seed1
  )

  high.embryo <- tessera::Embryo(
    n.cells = 200,
    n.chrs = 1,
    prop.aneuploid = p2,
    dispersal = d2,
    concordance = 1,
    rng.seed = seed2
  )

  low.biopsies <- tessera::takeAllBiopsies(low.embryo, biopsy.size = 5, chromosome = 1)
  high.biopsies <- tessera::takeAllBiopsies(high.embryo, biopsy.size = 5, chromosome = 1)

  combs <- expand.grid(low = low.biopsies, high = high.biopsies)
  combs$diff <- combs$low - combs$high
  calc.is.correct <- function(x) {
    if (x < 0) {
      return("Correct rank")
    }
    if (x == 0) {
      return("No rank")
    }
    return("Incorrect rank")
  }
  combs$rank <- sapply(combs$diff, calc.is.correct)

  # Return the total number of instances in each category
  c(
    "Correct rank" = length(combs$rank[combs$rank == "Correct rank"]),
    "Incorrect rank" = length(combs$rank[combs$rank == "Incorrect rank"]),
    "No rank" = length(combs$rank[combs$rank == "No rank"]),
    # Split the ties 50:50 between correct and incorrect
    "Adj.correct.rank" = (length(combs$rank[combs$rank == "Correct rank"]) + (length(combs$rank[combs$rank == "No rank"]) / 2)),
    "Adj.incorrect.rank" = (length(combs$rank[combs$rank == "Incorrect rank"]) + (length(combs$rank[combs$rank == "No rank"]) / 2))
  )
}


#' Calculate the ranks for a given combination of aneuploidy and dispersals
#' over 100 replicates
#'
#' @param p1 the fractional aneuploidy of embryo 1
#' @param p2 the fractional aneuploidy of embryo 2
#' @param d1 the dispersal of embryo 1
#' @param d2 the dispersal of embryo 2
#'
#' @return a named vector with the percentage of biopsies in each rank class over the 100 replicates
calculate.rank.combo <- function(p1, p2, d1, d2) {
  p <- mcmapply(calculate.ranks,
    p1 = p1, p2 = p2, d1 = d1, d2 = d2, seed1 = 0:100,
    seed2 = 0:100, mc.cores = N.CORES, SIMPLIFY = T
  )

  # Total number of biopsies ranked
  total <- sum(p["Correct rank", ]) + sum(p["Incorrect rank", ]) + sum(p["No rank", ])

  # Aggregated for combo
  c(
    "Low.aneu" = p1,
    "High.aneu" = p2,
    "Low.disp" = d1,
    "High.disp" = d2,
    "Correct.rank" = sum(p["Correct rank", ]) / total * 100,
    "Incorrect.rank" = sum(p["Incorrect rank", ]) / total * 100,
    "No.rank" = sum(p["No rank", ]) / total * 100,
    "Adj.correct.rank" = sum(p["Adj.correct.rank", ]) / total * 100,
    "Adj.incorrect.rank" = sum(p["Adj.incorrect.rank", ]) / total * 100
  )
}


# Calculate if results not available
out.file <- "data/Rank_results.csv"
if (!file.exists(out.file)) {
  aneuploidies <- expand.grid(
    p1 = seq(0, 1, 0.01),
    p2 = seq(0, 1, 0.01),
    d1 = seq(0, 1, 0.5),
    d2 = seq(0, 1, 0.5)
  )
  # Note - if allowing equal values, need to account for doubling of p1=p2 values
  aneuploidies <- aneuploidies[aneuploidies$p1 <= aneuploidies$p2, ]

  result <- mcmapply(calculate.rank.combo,
    p1 = aneuploidies$p1,
    p2 = aneuploidies$p2,
    d1 = aneuploidies$d1,
    d2 = aneuploidies$d2,
    mc.cores = N.CORES,
    SIMPLIFY = T
  )


  output <- as.data.frame(t(result))
  output$Aneu.diff <- round(output$High.aneu - output$Low.aneu, digits = 3)

  write.csv(output, file = "data/Rank_results.csv", quote = F, row.names = F, col.names = T)
}


# We also want to consider multiple embryos. How well can we select the best n
# embryos from a pool of m embryos?

# For i replicates at a given dispersal:
#   create a pool of m embryos with differing aneuploidies
#   biopsy each embryo and rank
#   select the best n embryos from biopsy
#   compare to the true rank order (how many matches to the true n?)
#     yields a percent correctly selected for the pool
# show the distribution of correct percents per dispersal
out.file <- "data/Rank_pool_results.csv"
if (!file.exists(out.file)) {
  POOL.SIZES <- seq(4, 10, 2) # number of embryos in the pool
  replicates <- 100 # number of pools to make for each dispersal

  combos <- expand.grid(
    "pool.size" = POOL.SIZES,
    "best.size" = 2:max(POOL.SIZES), # number of embryos to select for transfer
    "disp" = DISPERSAL.RANGE
  ) %>%
    dplyr::filter(best.size < pool.size)

  #' Calculate the percentage of embryos correctly selected from a pool via biopsy results.
  #' The pool contains embryos with randomly generated aneuploidies in each replicate
  #'
  #' @param pool.size number of embryos in the pool
  #' @param best.size number of embryos to select for transfer
  #' @param disp # dipsersal of the aneuploid cells in the embryo
  #'
  #' @return a data frame containing the percentage of correctly ranked embryos over 100 replicates
  calc.percents <- function(pool.size, best.size, disp) {
    cat("Pool:", pool.size, "  Best:", best.size, "  Disp:", disp, "\n")

    # Create an embryo with the given aneuploidy and seed, and take one random biopsy
    get.biopsy <- function(a, s) {
      e <- tessera::Embryo(n.cells = 200, n.chrs = 1, prop.aneuploid = a, dispersal = disp, rng.seed = s)
      tessera::takeBiopsy(e)
    }

    # find the percentage of embryos correctly selected by biopsy from a pool
    get.pct.correct <- function(i) {

      # different set of embryos for each iteration
      # biopsy will also be different even if embryo is reused
      seeds <- i:(i + (pool.size - 1))

      # Create a new random set of aneuploidies for each iteration
      aneuploidies <- runif(pool.size, min = 0, max = 1)
      biopsies <- mapply(get.biopsy, a = aneuploidies, s = seeds)

      # Rank the embryos by aneuploidy
      best.from.biopsy <- which(biopsies %in% sort(biopsies)[1:best.size])
      best.from.reality <- which(aneuploidies %in% sort(aneuploidies)[1:best.size])

      # Calculate the percent of embryos correctly selected for transfer by biopsy
      pct.corrrect <- sum(best.from.biopsy %in% best.from.reality) / best.size * 100
      pct.corrrect
    }

    correct.pcts <- sapply(1:replicates, get.pct.correct)

    data.frame(
      "disp" = disp,
      "pct.correct" = correct.pcts,
      "pool.size" = pool.size,
      "best.size" = best.size
    )
  }


  pool.results <- do.call(rbind, mcmapply(calc.percents,
    pool.size = combos$pool.size,
    best.size = combos$best.size,
    disp = combos$disp,
    SIMPLIFY = F,
    mc.cores = 20
  ))

  write.csv(pool.results, file = out.file, quote = F, row.names = F, col.names = T)
  # Figures generated in Figure_xxxx_rank_pool.R
}
