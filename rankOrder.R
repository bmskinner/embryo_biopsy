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
out.file <- "Rank_results.csv"
if (!file.exists(out.file)) {
  aneuploidies <- expand.grid(
    p1 = seq(0, 1, 0.05),
    p2 = seq(0, 1, 0.05),
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
  # head(output)

  write.csv(output, file = "Rank_results.csv", quote = F, row.names = F, col.names = T)
}

output <- read.csv("Rank_results.csv", header = T)

correct.plot <- ggplot(output, aes(x = Low.aneu, y = High.aneu, fill = Correct.rank)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0, 100)) +
  scale_y_continuous(
    sec.axis = sec_axis(~., name = "Embryo one dispersal", breaks = NULL, labels = NULL)
  ) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "Embryo two dispersal", breaks = NULL, labels = NULL)) +
  labs(x = "Embryo one aneuploidy", y = "Embryo two aneuploidy") +
  facet_grid(Low.disp ~ High.disp) +
  theme_classic()
save.double.width(correct.plot, "Figure xxxx - Rank_output_correct", 170)


no.rank.plot <- ggplot(output, aes(x = Low.aneu, y = High.aneu, fill = No.rank)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, 100)) +
  scale_y_continuous(
    sec.axis = sec_axis(~., name = "Embryo one dispersal", breaks = NULL, labels = NULL)
  ) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "Embryo two dispersal", breaks = NULL, labels = NULL)) +
  labs(x = "Embryo one aneuploidy", y = "Embryo two aneuploidy") +
  facet_grid(Low.disp ~ High.disp) +
  theme_classic()
save.double.width(no.rank.plot, "Figure xxxx - Rank_output_none", 170)

incorrect.rank.plot <- ggplot(output, aes(x = Low.aneu, y = High.aneu, fill = Incorrect.rank)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, 50)) +
  scale_y_continuous(
    sec.axis = sec_axis(~., name = "Embryo one dispersal", breaks = NULL, labels = NULL)
  ) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "Embryo two dispersal", breaks = NULL, labels = NULL)) +
  labs(x = "Embryo one aneuploidy", y = "Embryo two aneuploidy") +
  facet_grid(Low.disp ~ High.disp) +
  theme_classic()
save.double.width(incorrect.rank.plot, "Figure xxxx - Rank_output_incorrect", 170)


# Aggregate results by difference in aneuploidy and calculate means, SDs
filt <- output %>%
  # dplyr::filter(Low.disp == 0 & High.disp==0) %>%
  dplyr::group_by(Aneu.diff, Low.disp, High.disp) %>%
  dplyr::summarise(
    Mean.correct = mean(Correct.rank),
    SD.correct = sd(Correct.rank),
    Mean.incorrect = mean(Incorrect.rank),
    SD.incorrect = sd(Incorrect.rank),
    Mean.no.rank = mean(No.rank),
    SD.no.rank = sd(No.rank),
    Mean.adj.correct = mean(Adj.correct.rank),
    SD.adj.correct = sd(Adj.correct.rank),
    Mean.adj.incorrect = mean(Adj.incorrect.rank),
    SD.adj.incorrect = sd(Adj.incorrect.rank)
  )


# Show correct ranks
out.correct.plot <- ggplot(filt, aes(x = Aneu.diff, y = Mean.correct)) +
  annotate("rect", xmin = 0, xmax = 0.2, ymin = 0, ymax = Inf, fill = "lightgray") +
  geom_hline(yintercept = 50, col = "black") +
  geom_point() +
  geom_errorbar(aes(ymin = Mean.correct - SD.correct, ymax = Mean.correct + SD.correct), size = 0.5) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 20),
    sec.axis = sec_axis(~., name = "Embryo one dispersal", breaks = NULL, labels = NULL)
  ) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "Embryo two dispersal", breaks = NULL, labels = NULL)) +
  labs(x = "Difference in aneuploidy between embryos", y = "Percent correctly ranked") +
  facet_grid(Low.disp ~ High.disp) +
  theme_classic() +
  theme(
    axis.line.y = element_line(),
    panel.grid.major.y = element_line()
  )
save.double.width(out.correct.plot, "Figure xxxx - Rank_correct", 170)

# Show correct, incorrect and no ranks on single plot
out.combined.plot <- ggplot(filt, aes(x = Aneu.diff)) +
  annotate("rect", xmin = 0, xmax = 0.2, ymin = 0, ymax = Inf, fill = "lightgray") +
  geom_hline(yintercept = 50, col = "black") +
  geom_point(aes(y = Mean.correct), col = "blue") +
  geom_errorbar(aes(ymin = Mean.correct - SD.correct, ymax = Mean.correct + SD.correct), size = 0.5, col = "blue") +
  geom_point(aes(y = Mean.incorrect), col = "red", alpha = 0.5) +
  geom_errorbar(aes(ymin = Mean.incorrect - SD.incorrect, ymax = Mean.incorrect + SD.incorrect), size = 0.5, col = "red", alpha = 0.5) +
  geom_point(aes(y = Mean.no.rank), col = "black", alpha = 0.5) +
  geom_errorbar(aes(ymin = Mean.no.rank - SD.no.rank, ymax = Mean.no.rank + SD.no.rank), size = 0.5, col = "black", alpha = 0.5) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 20),
    sec.axis = sec_axis(~., name = "Embryo one dispersal", breaks = NULL, labels = NULL)
  ) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "Embryo two dispersal", breaks = NULL, labels = NULL)) +
  labs(x = "Difference in aneuploidy between embryos", y = "Percent correctly ranked") +
  facet_grid(Low.disp ~ High.disp) +
  theme_classic() +
  theme(
    axis.line.y = element_line(),
    panel.grid.major.y = element_line()
  )
save.double.width(out.combined.plot, "Figure xxxx - Ranks_all", 170)

# Split the values for no.rank between the correct and incorrect
# Can't get error bars direct from this - recalculate from original values
out.split.plot <- ggplot(filt, aes(x = Aneu.diff)) +
  geom_hline(yintercept = 50, col = "black") +
  geom_point(aes(y = Mean.adj.incorrect), col = "red") +
  geom_errorbar(aes(ymin = Mean.adj.incorrect - SD.adj.incorrect, ymax = Mean.adj.incorrect + SD.adj.incorrect), size = 0.5, col = "red") +
  geom_point(aes(y = Mean.adj.correct), col = "blue", alpha = 0.5) +
  geom_errorbar(aes(ymin = Mean.adj.correct - SD.adj.correct, ymax = Mean.adj.correct + SD.adj.correct), size = 0.5, col = "blue", alpha = 0.5) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 20),
    sec.axis = sec_axis(~., name = "Embryo one dispersal", breaks = NULL, labels = NULL)
  ) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "Embryo two dispersal", breaks = NULL, labels = NULL)) +
  labs(x = "Difference in aneuploidy between embryos", y = "Percent correctly ranked") +
  facet_grid(Low.disp ~ High.disp) +
  theme_classic() +
  theme(
    axis.line.y = element_line(),
    panel.grid.major.y = element_line()
  )
