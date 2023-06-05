# Figure S1 - 50% threshold for low-level / high-level aneuploidy

library(tessera)
library(tidyverse)
library(patchwork)
library(svglite)
library(data.table)

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

# Filter to 5-cell biopsy and calculate matching percentage
plot.values <- agg.values %>%
  dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size) %>%
  dplyr::summarise(
    pct_equal_match = mean(f_equal_match) * 100,
    sd_equal_match = sd(f_equal_match) * 100
  )

p <- ggplot(plot.values, aes(x = Aneuploidy * 100, y = Dispersal, fill = pct_equal_match)) +
  geom_raster() +
  geom_vline(xintercept = 19.5, col = "white") +
  geom_vline(xintercept = 39.5, col = "white") +
  geom_vline(xintercept = 80.5, col = "white") +
  labs(
    x = "Embryo aneuploidy",
    y = "Embryo dispersal",
    fill = "Percent of\nbiopsies in\ncorrect class"
  ) +
  coord_cartesian(ylim = c(0, 1.05)) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_fill_viridis_c(limits = c(0, 100)) +
  facet_wrap(~Biopsy_size, ncol = 4) +
  theme_classic() +
  theme(axis.title = element_text(size = 9))

save.double.width(p, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_S1_equal_classes"), height = 120)
