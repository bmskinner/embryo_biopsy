#  Figure 4 - like Fig 3 faceted over all biopsy sizes
library(tessera)
library(tidyverse)
library(patchwork)
library(svglite)
library(data.table)
library(parallel)

source("parameters.R")
source("functions.R")

# Read in the aggregate values
agg.values <- do.call(rbind, mclapply(list.files(
  path = AGGREGATE.DATA.PATH,
  pattern = "merged", full.names = T
),
fread,
header = T, mc.cores = N.CORES
))

# calculate matching percentage
plot.values <- agg.values %>%
  dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size) %>%
  dplyr::summarise(
    pct_equal_match = mean(f_equal_match) * 100,
    sd_equal_match = sd(f_equal_match) * 100
  )

p <- ggplot(plot.values, aes(x = Aneuploidy * 100, y = Dispersal, fill = pct_equal_match)) +
  geom_raster() +
  geom_vline(xintercept = 19.5, col = "white") +
  geom_vline(xintercept = 49.5, col = "white") +
  geom_vline(xintercept = 80.5, col = "white") +
  labs(
    x = "Embryo aneuploidy",
    y = "Embryo dispersal",
    fill = "Percent of\nbiopsies in\ncorrect class"
  ) +
  coord_cartesian(ylim = c(0, 1.05)) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_fill_viridis_c(limits = c(0, 100)) +
  facet_wrap(~Biopsy_size, ncol = 4, labeller = as_labeller(function(i) paste(i, "cells"))) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 9),
    panel.grid = element_blank()
  )

save.double.width(p, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_4_heatmap"), height = 120)
