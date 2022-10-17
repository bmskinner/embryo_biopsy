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
    "No rank" = length(combs$rank[combs$rank == "No rank"]) / nrow(combs) * 100,
    # Split the ties 50:50 between correct and incorrect
    "Adj.correct.rank" = (length(combs$rank[combs$rank == "Correct rank"])+(length(combs$rank[combs$rank == "No rank"])/2)) / nrow(combs) * 100,
    "Adj.incorrect.rank" = (length(combs$rank[combs$rank == "Incorrect rank"])+(length(combs$rank[combs$rank == "No rank"])/2)) / nrow(combs) * 100
  )
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
    "No rank" = mean(p["No rank", ]),
    "Adj.correct.rank" = mean(p["Adj.correct.rank", ]),
    "Adj.incorrect.rank" = mean(p["Adj.incorrect.rank", ])
  )
}


# Calculate if results not available
out.file = "Rank_results.csv"
if(!file.exists(out.file)){
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
  head(output)
  
  write.csv(output, file = "Rank_results.csv", quote = F, row.names = F, col.names = T)
}

output <- read.csv("Rank_results.csv", header = T)
output$Aneu.diff <- round(output$High.aneu - output$Low.aneu, digits = 3)


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


filt <- output %>%
  # dplyr::filter(Low.disp == 0 & High.disp==0) %>%
  dplyr::group_by(Aneu.diff, Low.disp, High.disp) %>%
  dplyr::summarise(
    Mean.correct = mean(Correct.rank),
    SD.correct = sd(Correct.rank),
    Mean.incorrect = mean(Incorrect.rank),
    SD.incorrect = sd(Incorrect.rank),
    Mean.no.rank = mean(No.rank),
    SD.no.rank = sd(No.rank)
  ) # %>%
# dplyr::select(Aneu.diff, Mean.correct, SD.correct) %>%
# dplyr::distinct()

# filt %>% tidyr::pivot_longer(cols=Mean.correct:SD.no.rank, names_to = "Type", values_to = "Value")

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

out.combined.plot <- ggplot(filt, aes(x = Aneu.diff)) +
  annotate("rect", xmin = 0, xmax = 0.2, ymin = 0, ymax = Inf, fill = "lightgray") +
  geom_hline(yintercept = 50, col = "black") +
  geom_point(aes(y = Mean.correct), col='blue') +
  geom_errorbar(aes(ymin = Mean.correct - SD.correct, ymax = Mean.correct + SD.correct), size = 0.5, col='blue') +
  geom_point(aes(y=Mean.incorrect), col='red', alpha=0.5)+
  geom_errorbar(aes(ymin = Mean.incorrect - SD.incorrect, ymax = Mean.incorrect + SD.incorrect), size = 0.5, col='red', alpha=0.5) +
  geom_point(aes(y=Mean.no.rank), col='black', alpha=0.5)+
  geom_errorbar(aes(ymin = Mean.no.rank - SD.no.rank, ymax = Mean.no.rank + SD.no.rank), size = 0.5, col='black', alpha=0.5) +
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
save.double.width(out.combined.plot, "Figure xxxx - Rank_combined", 170)

# Split the values for no.rank between the correct and incorrect
# Can't get error bars direct from this - recalculate from original values
out.split.plot <- ggplot(filt, aes(x = Aneu.diff)) +
  annotate("rect", xmin = 0, xmax = 0.2, ymin = 0, ymax = Inf, fill = "lightgray") +
  geom_hline(yintercept = 50, col = "black") +
  geom_point(aes(y = Mean.correct+(Mean.no.rank/2)), col='blue') +
  # geom_errorbar(aes(ymin = Mean.correct - SD.correct, ymax = Mean.correct + SD.correct), size = 0.5, col='blue') +
  geom_point(aes(y=Mean.incorrect+(Mean.no.rank/2)), col='red', alpha=0.5)+
  # geom_errorbar(aes(ymin = Mean.incorrect - SD.incorrect, ymax = Mean.incorrect + SD.incorrect), size = 0.5, col='red', alpha=0.5) +
  # geom_point(aes(y=Mean.no.rank), col='black', alpha=0.5)+
  # geom_errorbar(aes(ymin = Mean.no.rank - SD.no.rank, ymax = Mean.no.rank + SD.no.rank), size = 0.5, col='black', alpha=0.5) +
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
