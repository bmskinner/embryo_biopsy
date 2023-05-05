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
  dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size) %>%
  dplyr::summarise(
    MeanError = mean(diff_to_embryo) * 100,
    SDError = sd(diff_to_embryo) * 100
  )

fig2a <- ggplot(aneu.data, aes(x = Aneuploidy * 100, y = MeanError)) +
  geom_point() +
  geom_ribbon(aes(ymin = MeanError - SDError, ymax = MeanError + SDError), alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  labs(x = "Embryo aneuploidy", y = "Percent difference between\n5-cell biopsy and embryo") +
  theme_bw()

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
  dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size) %>%
  dplyr::summarise(
    MeanError = mean(diff_to_embryo),
    SDError = sd(diff_to_embryo)
  )

fig2b <- ggplot(disp.data, aes(x = Dispersal, y = MeanError)) +
  geom_point() +
  geom_ribbon(aes(ymin = MeanError - SDError, ymax = MeanError + SDError), alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  labs(x = "Embryo dispersal", y = "Percent difference between\n5-cell biopsy and embryo") +
  theme_bw()

fig2 <- fig2a + fig2b + plot_annotation(tag_levels = c("A"))
save.double.width(fig2, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_2"), height = 85)
