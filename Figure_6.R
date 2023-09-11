# Figure 6: Predictive heatmap and column plot for 5 cell biopsy

source("parameters.R")
source("functions.R")

library(png)
library(grid)
library(gridExtra)

# # Read the saved raw values for selected dispersals at 200 cell embryo
raw.values <- do.call(rbind, mclapply(list.files(
  path = RAW.DATA.PATH,
  pattern = "raw_values_e200_a.*_d(0|0.5|1).csv", full.names = T
),
fread,
header = T,
mc.cores = N.CORES
)) %>%
  dplyr::filter(Biopsy_size == 5)

fig6A.data <- raw.values %>%
  rowwise() %>%
  mutate(f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size * 100))) %>%
  select(-starts_with("V")) %>%
  tidyr::unnest(f_aneuploid) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(f_aneuploid, Aneuploidy, Dispersal) %>%
  dplyr::summarise(Count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(f_aneuploid, Dispersal) %>%
  dplyr::mutate(CountBiopsy = sum(Count)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(f_aneuploid, Aneuploidy, Dispersal) %>%
  dplyr::mutate(PctTotal = Count / CountBiopsy * 100)


# Blank canvas
zero.data <- expand.grid(
  Aneuploidy = ANEUPLOIDY.RANGE,
  f_aneuploid = seq(0, 100, 20),
  PctTotal = 0
)

# Part A - heatmap
fig6A.plot <- plot.equal.heatmap(fig6A.data, zero.data) + plot_annotation(tag_levels = list(c("A")))


# Calculate the proportion of embryos matching biopsy class
fig6B.data <- raw.values %>%
  rowwise() %>%
  mutate(n_aneuploid = list(as.double(c_across(starts_with("V"))))) %>%
  select(-starts_with("V")) %>%
  tidyr::unnest(n_aneuploid) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    EmbryoClass = to.equal.class(Aneuploidy),
    BiopsyClass = to.equal.class(n_aneuploid / Biopsy_size)
  ) %>%
  dplyr::group_by(Dispersal, n_aneuploid, EmbryoClass, BiopsyClass, Biopsy_size) %>%
  dplyr::summarise(Count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(n_aneuploid, Dispersal) %>%
  dplyr::mutate(
    Total = sum(Count),
    PctTotal = Count / Total * 100,
    IsCorrect = EmbryoClass == BiopsyClass
  ) %>%
  dplyr::filter(IsCorrect)

# Part B
fig6B.plot <- plot.columns(fig6B.data) + plot_annotation(tag_levels = list(c("B")))

save.single.width(fig6A.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_6A_predictive"), height = 85)
save.single.width(fig6B.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_6B_predictive"), height = 85)


################################################################################

# Column plots do not have a legend, and so should have plot boundaries that are
# higher than heatmaps. Patchwork would create too much whitespace.
# Rendered each panel to png, then combine.

fig6A.png <- readPNG(paste0(FIGURE.OUTPUT.DIR, "/Figure_6A_predictive.png"))
fig6B.png <- readPNG(paste0(FIGURE.OUTPUT.DIR, "/Figure_6B_predictive.png"))

# setup plot

png(filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_6_predictive.png"), width = 170, height = 85, units = "mm", res = 300)
par(mai = rep(0, 4)) # no margins

# do the plotting
plot(NA, xlim = 0:1, ylim = 0:1, bty = "n", axes = 0, xaxs = "i", yaxs = "i")
rasterImage(fig6A.png, 0, 0, 0.5, 1)
rasterImage(fig6B.png, 0.5, 0, 1, 1)
dev.off()



# Create all embryo size plots for supplement


# read.data <- function(x) {
#   fread(x, header = T) %>%
#     dplyr::filter(Biopsy_size == 5) %>%
#     rowwise() %>%
#     mutate(f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size * 100))) %>%
#     select(-starts_with("V")) %>%
#     tidyr::unnest(f_aneuploid)
# }
#
# all.files <- list.files(
#   path = RAW.DATA.PATH,
#   pattern = "raw_values_e.*_a.*_d(0|0.5|1).csv", full.names = T)
#
# all.values <- do.call(rbind, mclapply(all.files,
# fread,
# header = T,
# mc.cores = N.CORES
# ))
#
# fig.s1.data <- all.values %>%
#   rowwise() %>%
#   mutate(f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size * 100))) %>%
#   select(-starts_with("V")) %>%
#   tidyr::unnest(f_aneuploid) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(f_aneuploid, Aneuploidy, Dispersal, Embryo_size, Biopsy_size) %>%
#   dplyr::summarise(Count = dplyr::n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(f_aneuploid, Dispersal, Embryo_size, Biopsy_size) %>%
#   dplyr::mutate(CountBiopsy = sum(Count)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(f_aneuploid, Aneuploidy, Dispersal, Embryo_size, Biopsy_size) %>%
#   dplyr::mutate(PctTotal = Count / CountBiopsy * 100)
#
#
# calc.column.data <- function(data, class.function) {
#   data %>%
#     dplyr::mutate(
#       EmbryoClass = class.function(Aneuploidy),
#       BiopsyClass = class.function(n_aneuploid / Biopsy_size)
#     ) %>%
#     dplyr::group_by(Dispersal, n_aneuploid, Embryo_size, Biopsy_size, EmbryoClass, BiopsyClass) %>%
#     dplyr::summarise(Count = sum(Total_biopsies)) %>%
#     dplyr::ungroup() %>%
#     dplyr::group_by(Embryo_size, Biopsy_size, n_aneuploid, Dispersal) %>%
#     dplyr::mutate(
#       Total = sum(Count),
#       PctTotal = Count / Total * 100,
#       IsCorrect = EmbryoClass == BiopsyClass
#     ) %>%
#     dplyr::filter(IsCorrect)
# }
#
#
# make.column.plot <- function(data, biopsy.size) {
#   ggplot(data, aes(x = n_aneuploid / biopsy.size * 100, y = PctTotal, fill = IsCorrect)) +
#     geom_hline(yintercept = 50) +
#     geom_col(position = "stack") +
#     scale_fill_manual(values = c("dark green")) +
#     scale_x_continuous(
#       breaks = seq(0, 100, 20),
#       sec.axis = sec_axis(~., name = "Dispersal of aneuploid cells", breaks = NULL, labels = NULL)
#     ) +
#     scale_y_continuous(
#       breaks = seq(0, 100, 20),
#       sec.axis = sec_axis(~., name = "Embryo size", breaks = NULL, labels = NULL)
#     ) +
#     coord_cartesian(ylim = c(0, 100)) +
#     labs(
#       y = "Percentage of embryos\nmatching biopsy class",
#       x = "Biopsy aneuploidy"
#     ) +
#     facet_grid(Embryo_size ~ Dispersal) +
#     theme_bw() +
#     theme(legend.position = "none")
# }
#
# # Make two biopsy column plot
# col.data <- calc.column.data(plot.values, to.equal.class)
# col.plot <- make.two.biopsy.column.plot(col.data.equal, 10)
#
# save.double.width(col.plot.equal,
#                   filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_S1_predictive_columns_biopsy_size"),
#                   height = 150
# )
