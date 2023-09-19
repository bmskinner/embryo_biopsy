# Common functions
library(tessera)
library(tidyverse)
library(svglite)
library(data.table)
library(parallel)
library(patchwork)

# Return the equal split class for a given fractional aneuploidy
to.equal.class <- function(f.aneuploidy) {
  case_when(
    f.aneuploidy < 0.20 ~ "Euploid",
    f.aneuploidy < 0.50 ~ "Low level",
    f.aneuploidy <= 0.80 ~ "High level",
    f.aneuploidy <= 1 ~ "Aneuploid"
  )
}

# Save a plot with arbirary dimensions (in mm)
save.plot <- function(plot, filename, width, height) {
  ggsave(plot, filename = paste0(filename, ".png"), dpi = 300, units = "mm", width = width, height = height)
  ggsave(plot, filename = paste0(filename, ".svg"), dpi = 300, units = "mm", width = width, height = height)
}

# Save a plot as single column width (85mm)
save.single.width <- function(plot, filename, height) {
  save.plot(plot, filename, 85, height)
}

# Save a plot as double column width (170mm)
save.double.width <- function(plot, filename, height) {
  save.plot(plot, filename, 170, height)
}

# Make a basic heatmap with no rectangle annotations
plot.heatmap <- function(data, zero.data) {
  ggplot(data, aes(x = f_aneuploid, y = Aneuploidy * 100, fill = PctTotal)) +
    geom_raster(data = zero.data) +
    geom_raster() +
    scale_fill_viridis_c() +
    labs(
      x = "Biopsy aneuploidy (%)",
      y = "Embryo aneuploidy (%)",
      fill = "Percentage\nof biopsies"
    ) +
    scale_y_continuous(breaks = seq(0, 100, 20)) +
    scale_x_continuous(breaks = seq(0, 100, 20)) +
    theme_classic() +
    facet_wrap(~Dispersal, ncol = 3) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "top",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 9),
      legend.box.spacing = unit(1.5, "mm"),
      legend.key.height = unit(3, "mm"),
      legend.title.align = 1,
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
}


plot.equal.heatmap <- function(data, zero.data) {
  plot.heatmap(data, zero.data) +
    geom_rect(
      xmin = -10, xmax = 10, ymin = -2.5, ymax = 19.5,
      fill = NA, col = "white", size = 1
    ) +
    geom_rect(
      xmin = 10, xmax = 50, ymin = 19.5, ymax = 49.5,
      fill = NA, col = "white", size = 1
    ) +
    geom_rect(
      xmin = 50, xmax = 90, ymin = 49.5, ymax = 80.5,
      fill = NA, col = "white", size = 1
    ) +
    geom_rect(
      xmin = 90, xmax = 110, ymin = 80.5, ymax = 102.5,
      fill = NA, col = "white", size = 1
    )
}

# Plot sclassification accuracy column charts
plot.columns <- function(data) {
  ggplot(data, aes(x = n_aneuploid / Biopsy_size * 100, y = PctTotal, fill = IsCorrect)) +
    geom_hline(yintercept = 50) +
    geom_col(position = "stack") +
    scale_fill_manual(values = c(BIOPSY.COLUMN.RGB)) + # dark grey, hint of green
    scale_x_continuous(breaks = seq(0, 100, 20)) +
    scale_y_continuous(breaks = seq(0, 100, 20)) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(
      y = "Percentage of embryos\nmatching biopsy class",
      x = "Biopsy aneuploidy (%)"
    ) +
    facet_wrap(~Dispersal, ncol = 3) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
}


# plot - the heat map
draw.ten.cell.biopsy.classes <- function(plot) {
  plot +
    geom_rect(
      xmin = -5, xmax = 15, ymin = -0.025, ymax = 0.175,
      fill = NA, col = "white", size = 1
    ) +
    geom_rect(
      xmin = 15, xmax = 45, ymin = 0.175, ymax = 0.475,
      fill = NA, col = "white", size = 1
    ) +
    geom_rect(
      xmin = 45, xmax = 85, ymin = 0.475, ymax = 0.825,
      fill = NA, col = "white", size = 1
    ) +
    geom_rect(
      xmin = 85, xmax = 105, ymin = 0.825, ymax = 1.025,
      fill = NA, col = "white", size = 1
    )
}
