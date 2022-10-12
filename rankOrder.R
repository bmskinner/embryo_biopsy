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

  c(
    "Correct rank" = length(combs$rank[combs$rank == "Correct rank"]) / nrow(combs) * 100,
    "Incorrect rank" = length(combs$rank[combs$rank == "Incorrect rank"]) / nrow(combs) * 100,
    "No rank" = length(combs$rank[combs$rank == "No rank"]) / nrow(combs) * 100
  )
  # table(combs$rank) / sum(table(combs$rank))*100
}

calculate.rank.combo <- function(p1, p2, d1, d2) {
  p <- mcmapply(calculate.ranks,
    p1 = p1, p2 = p2, d1 = d1, d2 = d2, seed1 = 0:100,
    seed2 = 0:100, mc.cores = N.CORES, SIMPLIFY = T
  )

  c(
    "Low aneu" = p1,
    "High aneu" = p2,
    "Low disp" = d1,
    "High disp" = d2,
    "Correct rank" = mean(p["Correct rank", ]),
    "Incorrect rank" = mean(p["Incorrect rank", ]),
    "No rank" = mean(p["No rank", ])
  )
}


aneuploidies <- expand.grid(
  p1 = seq(0, 1, 0.05),
  p2 = seq(0, 1, 0.05),
  d1 = seq(0, 1, 0.5),
  d2 = seq(0, 1, 0.5)
)
aneuploidies <- aneuploidies[aneuploidies$p1 < aneuploidies$p2, ]

result <- mcmapply(calculate.rank.combo,
  p1 = aneuploidies$p1,
  p2 = aneuploidies$p2,
  d1 = aneuploidies$d1,
  d2 = aneuploidies$d2,
  mc.cores = N.CORES,
  SIMPLIFY = T
)


output <- as.data.frame(t(result))
head(output)

write.csv(output, file="Rank_results.csv", quote = F, row.names = F, col.names = T)

ggplot(output, aes(x = `Low aneu`, y = `High aneu`, fill = `Correct rank`)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0, 100)) +
  facet_grid(`Low disp` ~ `High disp`)
ggsave(filename = "Rank_output.png")
