# Figure 7:  predictive heatmaps and column plots with merged class thresholds

source("parameters.R")
source("functions.R")

# # Read the saved raw values for selected dispersals at 200 cell embryo
raw.values = do.call(rbind, mclapply(list.files(path = RAW.DATA.PATH,
                                                pattern = "raw_values_e200_a.*_d(0|0.5|1).csv", full.names = T),
                                     fread,
                                     header = T,
                                     mc.cores = N.CORES)) %>%
  dplyr::filter(Biopsy_size == 5)

heatmap.data = raw.values %>% 
  rowwise %>%
  mutate(f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size *100))) %>%
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
  dplyr::mutate(PctTotal = Count/CountBiopsy * 100)


# Blank canvas
zero.data = expand.grid(Aneuploidy = ANEUPLOIDY.RANGE,
                        f_aneuploid = seq(0, 100, 20),
                        PctTotal = 0)

# Part A and C - heatmap
hmap.plot  = plot.merge.heatmap(heatmap.data, zero.data)
save.single.width(hmap.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_7A_predictive"), height = 85)

hmap.plot2  = plot.merge2.heatmap(heatmap.data, zero.data)
save.single.width(hmap.plot2, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_7C_predictive"), height = 85)

# Part B and D
# Calculate the proportion of embryos matching biopsy class
col.data = raw.values %>% 
  rowwise %>%
  mutate(n_aneuploid = list(as.double(c_across(starts_with("V"))))) %>%
  select(-starts_with("V")) %>%
  tidyr::unnest(n_aneuploid) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(EmbryoClass = to.merged.class(Aneuploidy),
                BiopsyClass = to.merged.class(n_aneuploid/Biopsy_size)) %>%
  dplyr::group_by(Dispersal, n_aneuploid, EmbryoClass, BiopsyClass) %>%
  dplyr::summarise(Count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(n_aneuploid, Dispersal) %>%
  dplyr::mutate(Total = sum(Count),
                PctTotal = Count / Total * 100,
                IsCorrect = EmbryoClass == BiopsyClass) %>%
  dplyr::filter(IsCorrect)


col.plot = plot.columns(col.data)

save.single.width(col.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_7B_predictive"), height = 85)

col.data2 = raw.values %>% 
  rowwise %>%
  mutate(n_aneuploid = list(as.double(c_across(starts_with("V"))))) %>%
  select(-starts_with("V")) %>%
  tidyr::unnest(n_aneuploid) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(EmbryoClass = to.merged.class.2(Aneuploidy),
                BiopsyClass = to.merged.class.2(n_aneuploid/Biopsy_size)) %>%
  dplyr::group_by(Dispersal, n_aneuploid, EmbryoClass, BiopsyClass) %>%
  dplyr::summarise(Count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(n_aneuploid, Dispersal) %>%
  dplyr::mutate(Total = sum(Count),
                PctTotal = Count / Total * 100,
                IsCorrect = EmbryoClass == BiopsyClass) %>%
  dplyr::filter(IsCorrect)


col.plot2 = plot.columns(col.data2)

save.single.width(col.plot2, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_7D_predictive"), height = 85)
