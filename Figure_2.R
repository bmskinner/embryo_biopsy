# Figure 2: Plots for effect of changing aneuploidy and dispersal

# The goal of this figure is to show how aneuploidy and dispersal both
# contribute to the error between biopsy and embryo, but they behave in different
# ways

source("parameters.R")
source("functions.R")

# # Read the saved raw values for selected dispersals at 200 cell embryo
raw.aneu.values <- do.call(rbind, mclapply(list.files(
  path = RAW.DATA.PATH,
  pattern = "raw_values_e200_a.*_d1.csv",
  full.names = T
),
fread,
header = T,
mc.cores = N.CORES
))

# Aggregate the raw data
aneu.data <- raw.aneu.values %>%
  dplyr::filter(Biopsy_size == 5) %>%
  rowwise() %>%
  mutate(f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size))) %>%
  select(-starts_with("V")) %>%
  ungroup() %>%
  unnest(f_aneuploid) %>%
  dplyr::mutate(diff_to_embryo = round(abs(f_aneuploid - Aneuploidy), digits = 2)) %>%
  dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size, diff_to_embryo) %>%
  dplyr::summarise(Count = dplyr::n())

fig2a <- ggplot(aneu.data, aes(x = Aneuploidy, y = diff_to_embryo, fill = Count)) +
  geom_tile(width = 0.01, height = 0.01) +
  scale_fill_viridis_c() +
  labs(y = "Difference between embryo and biopsy", fill = "Number of\nbiopsies") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "#440154"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# Create plot for dispersal values

raw.disp.values <- do.call(rbind, mclapply(list.files(
  path = RAW.DATA.PATH,
  pattern = "raw_values_e200_a0.2_d.*.csv",
  full.names = T
),
fread,
header = T,
mc.cores = N.CORES
))

disp.data <- raw.disp.values %>%
  dplyr::filter(Biopsy_size == 5) %>%
  rowwise() %>%
  mutate(f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size))) %>%
  select(-starts_with("V")) %>%
  ungroup() %>%
  unnest(f_aneuploid) %>%
  dplyr::mutate(diff_to_embryo = round(abs(f_aneuploid - Aneuploidy), digits = 2)) %>%
  dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size, diff_to_embryo) %>%
  dplyr::summarise(Count = dplyr::n())


fig2b <- ggplot(disp.data, aes(x = Dispersal, y = diff_to_embryo, fill = Count)) +
  geom_tile(width = 0.01, height = 0.01) +
  scale_fill_viridis_c() +
  labs(y = "Difference between embryo and biopsy", fill = "Number of\nbiopsies") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "#440154"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

fig2 <- fig2a + fig2b + plot_annotation(tag_levels = c("A"))
save.double.width(fig2, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_2_heatmap"), height = 85)

################################################################################

# %>%
#   dplyr::bind_rows(., all.vals) %>% # add in all data
#   dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size, diff_to_embryo) %>%
#   dplyr::summarise(Count = sum(Count)) %>% # recalculate counts
#   dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size) %>%
#   dplyr::mutate(
#     TotalBiopsies = sum(Count),
#     PctBiopsies = Count / TotalBiopsies * 100
#   ) %>%
#   dplyr::distinct() %>%
#   dplyr::arrange(Biopsy_size, Dispersal, Aneuploidy, diff_to_embryo) %>%
#   dplyr::select(-TotalBiopsies, -Count) %>%
#   dplyr::distinct() %>%
#   dplyr::mutate(CumPct = cumsum(PctBiopsies))


# Given raw input data, calculate the difference between the biopsy and
# embryo aneuploidies and summarise into cumulative percentage
# calc.step.accuracy <- function(data, all.vals) {
#   data %>%
#     rowwise() %>%
#     mutate(f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size))) %>%
#     select(-starts_with("V")) %>%
#     ungroup() %>%
#     unnest(f_aneuploid) %>%
#     dplyr::mutate(diff_to_embryo = round(abs(f_aneuploid - Aneuploidy), digits = 2)) %>%
#     dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size, diff_to_embryo) %>%
#     dplyr::summarise(Count = dplyr::n()) %>%
#     dplyr::bind_rows(., all.vals) %>% # add in all data
#     dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size, diff_to_embryo) %>%
#     dplyr::summarise(Count = sum(Count)) %>% # recalculate counts
#     dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size) %>%
#     dplyr::mutate(
#       TotalBiopsies = sum(Count),
#       PctBiopsies = Count / TotalBiopsies * 100
#     ) %>%
#     dplyr::distinct() %>%
#     dplyr::arrange(Biopsy_size, Dispersal, Aneuploidy, diff_to_embryo) %>%
#     dplyr::select(-TotalBiopsies, -Count) %>%
#     dplyr::distinct() %>%
#     dplyr::mutate(CumPct = cumsum(PctBiopsies))
# }
#
#
# # Filter down so the plot is not overcrowded
#
#
#
#
#
# step.aneuploidy <- calc.step.accuracy(raw.values[raw.values$Dispersal == 1 & Biopsy_size == 5, ], all.vals)
#
# # Calculate step values for dispersal
# # All values for filling out the steps
# all.vals <- expand.grid(
#   "Aneuploidy" = 0.2,
#   "Dispersal" = DISPERSAL.RANGE,
#   "Biopsy_size" = 5,
#   "diff_to_embryo" = round(seq(0, 1, 0.01), digits = 2),
#   "Count" = 0
# )
#
# # Read in the files needed for dispersal step calculation
# step.values <- do.call(rbind, mclapply(list.files(
#   path = RAW.DATA.PATH,
#   pattern = "raw_values_a0.2_d.*.csv", full.names = T
# ),
# fread,
# header = T,
# mc.cores = 1
# ))
#
# step.dispersal <- calc.step.accuracy(step.values[step.values$Aneuploidy == 0.2 & Biopsy_size == 5, ], all.vals)
#
#
# step.aneuploidy.filt <- step.aneuploidy[step.aneuploidy$Aneuploidy %in% seq(0, 0.5, 0.05), ]
#
# aneuploidy.step.plot <- ggplot(step.aneuploidy.filt, aes(
#   x = diff_to_embryo * 100, y = CumPct,
#   col = Aneuploidy * 100, group = Aneuploidy
# )) +
#   geom_step(size = 2) +
#   labs(
#     y = "Percent of biopsies",
#     x = "Percent difference\nbetween biopsy and embryo",
#     col = "Aneuploidy"
#   ) +
#   scale_colour_viridis_c() +
#   scale_x_continuous(breaks = seq(0, 100, 20)) +
#   scale_y_continuous(breaks = seq(0, 100, 20)) +
#   theme_classic() +
#   theme(
#     legend.position = c(0.85, 0.4),
#     legend.title = element_text(size = 9),
#     axis.title = element_text(size = 9)
#   )
#
# # Dispersal step plot
#
# dispersal.step.plot <- ggplot(step.dispersal, aes(
#   x = diff_to_embryo * 100, y = CumPct,
#   col = Dispersal, group = Dispersal
# )) +
#   geom_step(size = 2) +
#   labs(
#     y = "Percent of biopsies",
#     x = "Percent difference\nbetween biopsy and embryo"
#   ) +
#   scale_colour_viridis_c() +
#   scale_x_continuous(breaks = seq(0, 100, 20)) +
#   scale_y_continuous(breaks = seq(0, 100, 20)) +
#   theme_classic() +
#   theme(
#     legend.position = c(0.85, 0.4),
#     legend.title = element_text(size = 9),
#     axis.title = element_text(size = 9)
#   )
#
# fig2 <- aneuploidy.step.plot + dispersal.step.plot + plot_annotation(tag_levels = c("A"))
# save.double.width(fig2, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_2_step"), height = 85)
