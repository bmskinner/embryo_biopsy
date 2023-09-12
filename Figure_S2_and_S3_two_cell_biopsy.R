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

make.two.biopsy.column.plot <- function(data, biopsy.size) {
  ggplot(data, aes(x = n_aneuploid / biopsy.size * 100, y = PctTotal, fill = IsCorrect)) +
    geom_hline(yintercept = 50) +
    geom_col(position = "stack") +
    scale_fill_manual(values = c(BIOPSY.COLUMN.RGB)) +
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


b <- 5
filt.data <- in.data %>%
  dplyr::filter(Biopsy_size == b)

# Blank canvas
zero.data <- expand.grid(
  Aneuploidy = ANEUPLOIDY.RANGE,
  f_aneuploid = seq(0, 100, 100 / (b * 2)),
  PctTotal = 0
)

# Plot two biopsy heatmap
hmap.plot.equal <- make.two.biopsy.heatmap(filt.data, zero.data)
hmap.plot.equal <- draw.ten.cell.biopsy.classes(hmap.plot.equal)

save.double.width(hmap.plot.equal,
  filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_S2b_predictive_heatmap_two_biopsy_", b),
  height = 150
)

# Make two biopsy column plot
col.data.equal <- filt.data %>%
  dplyr::filter(Biopsy_size == b) %>%
  dplyr::mutate(
    EmbryoClass = to.equal.class(Aneuploidy),
    BiopsyClass = to.equal.class(n_aneuploid / Biopsy_size)
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

# Save column values to be used in two biopsy comparisons
write.csv(col.data.equal,
  file = paste0("data/2x5_cell_biopsy_predictive_columns.csv"),
  quote = F, row.names = F, col.names = T
)

col.plot.equal <- make.two.biopsy.column.plot(col.data.equal, b)

save.double.width(col.plot.equal,
  filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_S2_predictive_columns_two_biopsy_", b),
  height = 150
)

# how does this compare to a 10-cell or five cell biopsy?

data.1x10 <- fread("data/1x10_cell_biopsy_predictive_columns.csv", header=T)
data.2x5 <-  fread("data/2x5_cell_biopsy_predictive_columns.csv", header=T)

data.1x10$Type = "1x10"
data.2x5$Type = "2x5"
data.2x5$n_aneuploid = data.2x5$n_aneuploid*2

data.all <- rbind(data.1x10, data.2x5) %>%
  tidyr::pivot_wider(id_cols = c(n_aneuploid, Dispersal, EmbryoClass, BiopsyClass, Embryo_size), 
                     values_from = PctTotal, 
                     names_from = Type) %>%
  dplyr::mutate(Diff = `2x5`- `1x10` )


plot.2x5.vs.1x10 <- ggplot(data.all, aes(x = n_aneuploid* 10, y = Diff)) +
  geom_hline(yintercept = 0) +
  geom_col(position = "stack", fill = BIOPSY.COLUMN.RGB) +
  # scale_fill_manual(values = c(BIOPSY.COLUMN.RGB)) +
  scale_x_continuous(
    breaks = seq(0, 100, 20),
    sec.axis = sec_axis(~., name = "Dispersal of aneuploid cells", breaks = NULL, labels = NULL)
  ) +
  scale_y_continuous(
    breaks = seq(-20, 20, 5),
    sec.axis = sec_axis(~., name = "Embryo size", breaks = NULL, labels = NULL)
  ) +
  coord_cartesian(ylim = c(-15, 15)) +
  labs(
    y = "Accuracy of two 5-cell biopsies compared to one 10-cell biopsy (%)",
    x = "Biopsy aneuploidy (%)"
  ) +
  facet_grid(Embryo_size ~ Dispersal) +
  theme_bw() +
  theme(legend.position = "none")

save.double.width(plot.2x5.vs.1x10,
                  filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_S3_two_biopsy_vs_one"),
                  height = 170
)
