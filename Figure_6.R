# Figure 6: Predictive heatmap and column plot for 5 cell biopsy

source("parameters.R")
source("functions.R")

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
fig6A.plot <- plot.pgdis.heatmap(fig6A.data, zero.data)


# Calculate the proportion of embryos matching biopsy class
fig6B.data <- raw.values %>%
  rowwise() %>%
  mutate(n_aneuploid = list(as.double(c_across(starts_with("V"))))) %>%
  select(-starts_with("V")) %>%
  tidyr::unnest(n_aneuploid) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    EmbryoClass = to.pgdis.class(Aneuploidy),
    BiopsyClass = to.pgdis.class(n_aneuploid / Biopsy_size)
  ) %>%
  dplyr::group_by(Dispersal, n_aneuploid, EmbryoClass, BiopsyClass) %>%
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
fig6B.plot <- plot.columns(fig6B.data)

save.single.width(fig6A.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_6A_predictive"), height = 85)
save.single.width(fig6B.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_6B_predictive"), height = 85)


# Part C - precictive power with merged mosaic class

fig6c.data <- raw.values %>%
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

fig6C.plot <- plot.merge.heatmap(fig6c.data, zero.data)
save.single.width(fig6C.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_6C_predictive"), height = 85)

# Part D - column
# Calculate the proportion of embryos matching biopsy class
fig6D.data <- raw.values %>%
  rowwise() %>%
  mutate(n_aneuploid = list(as.double(c_across(starts_with("V"))))) %>%
  select(-starts_with("V")) %>%
  tidyr::unnest(n_aneuploid) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    EmbryoClass = to.merged.class(Aneuploidy),
    BiopsyClass = to.merged.class(n_aneuploid / Biopsy_size)
  ) %>%
  dplyr::group_by(Dispersal, n_aneuploid, EmbryoClass, BiopsyClass) %>%
  dplyr::summarise(Count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(n_aneuploid, Dispersal) %>%
  dplyr::mutate(
    Total = sum(Count),
    PctTotal = Count / Total * 100,
    IsCorrect = EmbryoClass == BiopsyClass
  ) %>%
  dplyr::filter(IsCorrect)


fig6D.plot <- plot.columns(fig6D.data)
save.single.width(fig6D.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_6D_predictive"), height = 85)
