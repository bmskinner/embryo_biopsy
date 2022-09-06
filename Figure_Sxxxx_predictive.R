# Figure Sxxxx - Predictive heatmap across embryo sizes, 5 cell biopsy
source("parameters.R")
source("functions.R")

# # Read the saved raw values for selected dispersals
raw.values = do.call(rbind, mclapply(list.files(path = RAW.DATA.PATH,
                                                pattern = "raw_values_e.*_a.*_d(0|0.5|1).csv", full.names = T),
                                     function(x){fread(x) %>% rowwise %>%
                                         mutate(f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size *100))) %>%
                                         select(-starts_with("V")) %>%
                                         tidyr::unnest(f_aneuploid) },
                                     header = T,
                                     mc.cores = N.CORES)) %>%
  dplyr::filter(Biopsy_size == 5)

heatmap.data = raw.values %>% 
  # rowwise %>%
  # mutate(f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size *100))) %>%
  # select(-starts_with("V")) %>%
  # tidyr::unnest(f_aneuploid) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(f_aneuploid, Aneuploidy, Dispersal, Embryo_size) %>%
  dplyr::summarise(Count = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(f_aneuploid, Dispersal) %>%
  dplyr::mutate(CountBiopsy = sum(Count)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(f_aneuploid, Aneuploidy, Dispersal, Embryo_size) %>%
  dplyr::mutate(PctTotal = Count/CountBiopsy * 100)


# Blank canvas
zero.data = expand.grid(Aneuploidy = ANEUPLOIDY.RANGE,
                        f_aneuploid = seq(0, 100, 20),
                        Embryo_size = EMBRYO.SIZES,
                        PctTotal = 0)


hmap.plot  = ggplot(heatmap.data, aes(x = f_aneuploid, y = Aneuploidy, fill=PctTotal))+
  geom_raster(data = zero.data)+
  geom_raster() +
  geom_rect(xmin=-10, xmax=10, ymin=-0.025, ymax=0.195, 
            fill=NA, col="white", size=1)+
  geom_rect(xmin=10, xmax=30, ymin=0.195, ymax=0.395, 
            fill=NA, col="white", size=1)+
  geom_rect(xmin=30, xmax=90, ymin=0.395, ymax=0.805, 
            fill=NA, col="white", size=1)+
  geom_rect(xmin=90, xmax=110, ymin=0.805, ymax=1.025, 
            fill=NA, col="white", size=1)+
  scale_fill_viridis_c() +
  labs(x = "Biopsy aneuploidy",
       y = "Embryo aneuploidy",
       fill = "Percentage\nof biopsies") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  theme_classic() + 
  facet_wrap(Embryo_size~Dispersal, ncol = 3)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

save.double.width(hmap.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_Sxxxx_predictive"), height = 120)