# Figure S1 and S2 - Predictive heatmap and columns across embryo sizes, 10 cell biopsy
source("parameters.R")
source("functions.R")

read.data <- function(x) {
  fread(x, header = T) %>%
    rowwise() %>%
    mutate(f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size * 100))) %>%
    select(-starts_with("V")) %>%
    tidyr::unnest(f_aneuploid)
}

# # Read the saved raw values for selected dispersals
raw.values <- do.call(rbind, mclapply(list.files(
  path = RAW.DATA.PATH,
  pattern = "raw_values_e.*_a.*_d(0|0.5|1).csv", full.names = T
),
read.data,
mc.cores = N.CORES
)) %>%
  dplyr::filter(Biopsy_size == 10)

heatmap.data <- raw.values %>%
  dplyr::ungroup() %>%
  dplyr::group_by(f_aneuploid, Aneuploidy, Dispersal, Embryo_size, Biopsy_size) %>%
  dplyr::summarise(Count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(f_aneuploid, Dispersal, Embryo_size, Biopsy_size) %>%
  dplyr::mutate(CountBiopsy = sum(Count)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(f_aneuploid, Aneuploidy, Dispersal, Embryo_size, Biopsy_size) %>%
  dplyr::mutate(PctTotal = Count / CountBiopsy * 100)

# Predictive heatmap
hmap.plot <- ggplot(heatmap.data, aes(x = f_aneuploid, y = Aneuploidy, fill = PctTotal)) +

  # Zero data background fill
  geom_rect(
    xmin = -5, xmax = 105, ymin = 0, ymax = 1,
    fill = VIRIDIS.MIN.RGB
  ) +
  geom_raster() +
  scale_fill_viridis_c() +
  labs(
    x = "Biopsy aneuploidy",
    y = "Embryo aneuploidy",
    fill = "Percentage\nof biopsies"
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    sec.axis = sec_axis(~., name = "Embryo size", breaks = NULL, labels = NULL)
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, 20),
    sec.axis = sec_axis(~., name = "Dispersal of aneuploid cells", breaks = NULL, labels = NULL)
  ) +
  theme_classic() +
  facet_grid(Embryo_size ~ Dispersal)

hmap.plot <- draw.ten.cell.biopsy.classes(hmap.plot)

save.double.width(hmap.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_S1b_10-cell_predictive_heatmap"), height = 170)

################################################################################

# Read the raw data for column calculations

################################################################################

read.class.data <- function(x) {
  fread(x, header = T) %>%
    dplyr::filter(Biopsy_size == 10) %>%
    rowwise() %>%
    mutate(n_aneuploid = list(as.double(c_across(starts_with("V"))))) %>%
    select(-starts_with("V")) %>%
    tidyr::unnest(n_aneuploid)
}

raw.class.data <- do.call(rbind, mclapply(list.files(
  path = RAW.DATA.PATH,
  pattern = "raw_values_e.*_a.*_d(0|0.5|1).csv", full.names = T
),
read.class.data,
mc.cores = N.CORES
))


# Correct class columns
col.data <- raw.class.data %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    EmbryoClass = to.equal.class(Aneuploidy),
    BiopsyClass = to.equal.class(n_aneuploid / Biopsy_size)
  ) %>%
  dplyr::group_by(Dispersal, n_aneuploid, EmbryoClass, BiopsyClass, Embryo_size) %>%
  dplyr::summarise(Count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(n_aneuploid, Dispersal, Embryo_size) %>%
  dplyr::mutate(
    Total = sum(Count),
    PctTotal = Count / Total * 100,
    IsCorrect = EmbryoClass == BiopsyClass
  ) %>%
  dplyr::filter(IsCorrect)

# Save column values to be used in two biopsy comparisons
write.csv(col.data,
  file = paste0("data/1x10_cell_biopsy_predictive_columns.csv"),
  quote = F, row.names = F, col.names = T
)


col.plot <- ggplot(col.data, aes(x = n_aneuploid * 10, y = PctTotal, fill = IsCorrect)) +
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

save.double.width(col.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_S1_10-cell_predictive_columns"), height = 170)
