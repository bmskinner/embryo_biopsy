# Figure 3: 5-cell biopsy heatmap
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
  dplyr::filter(Biopsy_size == 5) %>%
  dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size) %>%
  dplyr::summarise(
    pct_pgdis_match = mean(f_pgdis_match) * 100,
    sd_pgdis_match = sd(f_pgdis_match) * 100
  )

p <- ggplot(plot.values, aes(x = Aneuploidy * 100, y = Dispersal, fill = pct_pgdis_match)) +
  geom_raster() +
  geom_vline(xintercept = 19.5, col = "white") +
  geom_vline(xintercept = 39.5, col = "white") +
  geom_vline(xintercept = 80.5, col = "white") +
  annotate(geom = "text", label = "Euploid", x = 10, y = 1.06, size = 3) +
  annotate(geom = "text", label = "Low", x = 30, y = 1.06, size = 3) +
  annotate(geom = "text", label = "High", x = 60, y = 1.06, size = 3) +
  annotate(geom = "text", label = "Aneuploid", x = 91, y = 1.06, size = 3) +
  labs(
    x = "Embryo aneuploidy",
    y = "Embryo dispersal",
    fill = "Percent of biopsies\nin correct class"
  ) +
  coord_cartesian(ylim = c(0, 1.05)) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_fill_viridis_c(limits = c(0, 100)) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9),
    legend.box.spacing = unit(1.5, "mm"),
    legend.key.height = unit(3, "mm"),
    legend.title.align = 1,
    axis.title = element_text(size = 9),
    panel.grid = element_blank()
  )

save.single.width(p, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_3_heatmap"), height = 85)
