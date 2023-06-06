# Create two-biopsy figures
library(tessera)
library(parallel)
library(tidyverse)
library(data.table)
library(fs)

source("parameters.R")
source("functions.R")

# Read and analyse the saved data
in.files <- list.files(
  path = TWO.BIOPSY.DATA.PATH,
  pattern = paste0("two_biopsy_values_a.*_d(0|0.5|1).csv"), full.names = T
)

in.data <- do.call(rbind, mclapply(in.files,
  fread,
  header = T,
  mc.cores = N.CORES
))

in.data <- in.data %>%
  dplyr::group_by(Embryo_size, Dispersal, Biopsy_size) %>%
  dplyr::mutate(
    f_aneuploid = n_aneuploid / Biopsy_size * 100,
    Biopsies_from_combo = sum(Total_biopsies),
    PctTotal = Total_biopsies / Biopsies_from_combo * 100
  )

make.two.biopsy.heatmap <- function(data, zero.data) {
  ggplot(data, aes(x = f_aneuploid, y = Aneuploidy, fill = PctTotal)) +
    geom_raster(data = zero.data) +
    geom_raster() +
    scale_fill_viridis_c() +
    labs(
      x = "Biopsy aneuploidy",
      y = "Embryo aneuploidy",
      fill = "Percentage\nof biopsies"
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, 0.1),
      sec.axis = sec_axis(~., name = "Embryo size", breaks = NULL, labels = NULL)
    ) +
    scale_x_continuous(
      breaks = seq(0, 100, 20),
      sec.axis = sec_axis(~., name = "Dispersal of aneuploid cells", breaks = NULL, labels = NULL)
    ) +
    facet_grid(Embryo_size ~ Dispersal) +
    theme_bw() +
    theme(panel.grid = element_blank())
}

# plot - the heat map
# biopsy.size - the size of biopsy to annotate
draw.biopsy.pgdis.classes <- function(plot, biopsy.size) {
  plot +
    geom_rect(
      xmin = -5, xmax = 15, ymin = -0.025, ymax = 0.175,
      fill = NA, col = "white", size = 1
    ) +
    geom_rect(
      xmin = 15, xmax = 35, ymin = 0.175, ymax = 0.375,
      fill = NA, col = "white", size = 1
    ) +
    geom_rect(
      xmin = 35, xmax = 85, ymin = 0.375, ymax = 0.825,
      fill = NA, col = "white", size = 1
    ) +
    geom_rect(
      xmin = 85, xmax = 105, ymin = 0.825, ymax = 1.025,
      fill = NA, col = "white", size = 1
    )
}

draw.biopsy.merge.classes <- function(plot, biopsy.size) {
  plot +
    geom_rect(
      xmin = -5, xmax = 15, ymin = -0.025, ymax = 0.175,
      fill = NA, col = "white", size = 1
    ) +
    geom_rect(
      xmin = 15, xmax = 85, ymin = 0.175, ymax = 0.825,
      fill = NA, col = "white", size = 1
    ) +
    geom_rect(
      xmin = 85, xmax = 105, ymin = 0.825, ymax = 1.025,
      fill = NA, col = "white", size = 1
    )
}


calc.column.data <- function(data, class.function) {
  data %>%
    dplyr::filter(Biopsy_size == b) %>%
    dplyr::mutate(
      EmbryoClass = class.function(Aneuploidy),
      BiopsyClass = class.function(n_aneuploid / Biopsy_size)
    ) %>%
    dplyr::group_by(Dispersal, n_aneuploid, Embryo_size, EmbryoClass, BiopsyClass) %>%
    dplyr::summarise(Count = sum(Total_biopsies)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Embryo_size, n_aneuploid, Dispersal) %>%
    dplyr::mutate(
      Total = sum(Count),
      PctTotal = Count / Total * 100,
      IsCorrect = EmbryoClass == BiopsyClass
    ) %>%
    dplyr::filter(IsCorrect)
}

make.two.biopsy.column.plot <- function(data, biopsy.size) {
  ggplot(data, aes(x = n_aneuploid / biopsy.size * 100, y = PctTotal, fill = IsCorrect)) +
    geom_hline(yintercept = 50) +
    geom_col(position = "stack") +
    scale_fill_manual(values = c("dark green")) +
    scale_x_continuous(
      breaks = seq(0, 100, 20),
      sec.axis = sec_axis(~., name = "Dispersal of aneuploid cells", breaks = NULL, labels = NULL)
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, 20),
      sec.axis = sec_axis(~., name = "Embryo size", breaks = NULL, labels = NULL)
    ) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(
      y = "Percentage of embryos\nmatching biopsy class",
      x = "Biopsy aneuploidy"
    ) +
    facet_grid(Embryo_size ~ Dispersal) +
    theme_bw() +
    theme(legend.position = "none")
}


# for (b in BIOPSY.SIZES) {
b <- 5
filt.data <- in.data %>%
  dplyr::filter(Biopsy_size == b)

# Blank canvas
zero.data <- expand.grid(
  Aneuploidy = ANEUPLOIDY.RANGE,
  f_aneuploid = seq(0, 100, 100 / (b * 2)),
  PctTotal = 0
)

# Plot with PGDIS classes
hmap.plot.pgdis <- make.two.biopsy.heatmap(filt.data, zero.data)
hmap.plot.pgdis <- draw.biopsy.pgdis.classes(hmap.plot.pgdis, b)

save.double.width(hmap.plot.pgdis,
  filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_S2_predictive_heatmap_two_biopsy_", b),
  height = 150
)

# Plot with merge classes
hmap.plot.merge <- make.two.biopsy.heatmap(filt.data, zero.data)
hmap.plot.merge <- draw.biopsy.merge.classes(hmap.plot.merge, b)

save.double.width(hmap.plot.merge,
  filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_heatmap_two_biopsy_merge_", b),
  height = 150
)

# Make columns
col.data.pgdis <- calc.column.data(filt.data, to.pgdis.class)
col.plot.pgdis <- make.two.biopsy.column.plot(col.data.pgdis, b)

save.double.width(col.plot.pgdis,
  filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_S3_predictive_columns_two_biopsy_", b),
  height = 150
)

col.data.merge <- calc.column.data(filt.data, to.merged.class)
col.plot.merge <- make.two.biopsy.column.plot(col.data.merge, b)

save.double.width(col.plot.merge,
  filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_columns_two_biopsy_merge_", b),
  height = 150
)
