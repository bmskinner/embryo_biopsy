# Create two-biopsy figures
library(tessera)
library(parallel)
library(tidyverse)
library(data.table)
library(fs)

source("parameters.R")
source("functions.R")

# Read and analyse the saved data
in.files = list.files(path = TWO.BIOPSY.DATA.PATH,
                      pattern = paste0("two_biopsy"), full.names = T)

in.data = do.call(rbind, mclapply(in.files,
                                  fread,
                                  header = T,
                                  mc.cores = N.CORES))

for(b in BIOPSY.SIZES){
  
  filt.data = in.data %>% 
    dplyr::group_by(Embryo_size, Dispersal, Aneuploidy, Biopsy_size) %>%
    dplyr::filter(Biopsy_size == b) %>%
    dplyr::mutate(f_aneuploid         = n_aneuploid/Biopsy_size*100,
                  Biopsies_from_combo = sum(Total_biopsies),
                  PctTotal = Total_biopsies/Biopsies_from_combo*100)

  # Blank canvas
  zero.data = expand.grid(Aneuploidy = ANEUPLOIDY.RANGE,
                          f_aneuploid = seq(0, 100, 100/(b*2)),
                          PctTotal = 0)
  
  # Plot with PGDIS classes
  hmap.plot.pgdis = ggplot(filt.data, aes(x = f_aneuploid, y = Aneuploidy, fill=PctTotal))+
    geom_raster(data = zero.data)+
    geom_raster() +
    geom_rect(xmin=-5, xmax=15, ymin=-0.025, ymax=0.175, 
              fill=NA, col="white", size=1)+
    geom_rect(xmin=15, xmax=35, ymin=0.175, ymax=0.375, 
              fill=NA, col="white", size=1)+
    geom_rect(xmin=35, xmax=85, ymin=0.375, ymax=0.825, 
              fill=NA, col="white", size=1)+
    geom_rect(xmin=85, xmax=105, ymin=0.825, ymax=1.025, 
              fill=NA, col="white", size=1)+
    scale_fill_viridis_c() +
    labs(x = "Biopsy aneuploidy",
         y = "Embryo aneuploidy",
         fill = "Percentage\nof biopsies") +
    scale_y_continuous(breaks = seq(0, 1, 0.1),
                       sec.axis = sec_axis(~ . , name = "Embryo size",breaks = NULL, labels = NULL)) +
    scale_x_continuous(breaks = seq(0, 100, 20), 
                       sec.axis = sec_axis(~ . , name = "Dispersal of aneuploid cells",breaks = NULL, labels = NULL)) +
    facet_grid(Embryo_size~Dispersal)+
    theme_classic()

  save.double.width(hmap.plot.pgdis, 
    filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_heatmap_two_biopsy_pgdis_", b), 
    height = 150)
  
  # Plot with merge classes
  hmap.plot.merge = ggplot(filt.data, aes(x = f_aneuploid, y = Aneuploidy, fill=PctTotal))+
    geom_raster(data = zero.data)+
    geom_raster() +
    geom_rect(xmin=-5, xmax=15, ymin=-0.025, ymax=0.175, 
              fill=NA, col="white", size=1)+
    geom_rect(xmin=15, xmax=85, ymin=0.175, ymax=0.825, 
              fill=NA, col="white", size=1)+
    geom_rect(xmin=85, xmax=105, ymin=0.825, ymax=1.025, 
              fill=NA, col="white", size=1)+
    scale_fill_viridis_c() +
    labs(x = "Biopsy aneuploidy",
         y = "Embryo aneuploidy",
         fill = "Percentage\nof biopsies") +
    scale_y_continuous(breaks = seq(0, 1, 0.1),
                       sec.axis = sec_axis(~ . , name = "Embryo size",breaks = NULL, labels = NULL)) +
    scale_x_continuous(breaks = seq(0, 100, 20), 
                       sec.axis = sec_axis(~ . , name = "Dispersal of aneuploid cells",breaks = NULL, labels = NULL)) +
    facet_grid(Embryo_size~Dispersal)+
    theme_classic()

  save.double.width(hmap.plot.merge, 
    filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_heatmap_two_biopsy_merge_", b), 
    height = 150)

  col.data.pgdis = in.data %>% 
    dplyr::filter(Biopsy_size == b) %>%
    dplyr::mutate(EmbryoPGDISClass = to.pgdis.class(Aneuploidy),
                  BiopsyPGDISClass = to.pgdis.class(n_aneuploid/Biopsy_size)) %>%
    dplyr::group_by(Dispersal, n_aneuploid, Embryo_size, EmbryoPGDISClass, BiopsyPGDISClass) %>%
    dplyr::summarise(Count = sum(Total_biopsies)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Embryo_size, n_aneuploid, Dispersal) %>%
    dplyr::mutate(Total = sum(Count),
                  PctTotal = Count / Total * 100,
                  IsPGDISCorrect = EmbryoPGDISClass == BiopsyPGDISClass) %>%
    dplyr::filter(IsPGDISCorrect)
  
 col.plot.pgdis = ggplot(col.data.pgdis, aes(x = n_aneuploid/b*100, y= PctTotal,  fill=IsPGDISCorrect))+
    geom_hline(yintercept = 50) +
    geom_col(position="stack") +
    scale_fill_manual(values = c("dark green")) +
    scale_x_continuous(breaks = seq(0, 100, 20), 
                       sec.axis = sec_axis(~ . , name = "Dispersal of aneuploid cells",breaks = NULL, labels = NULL)) +
    scale_y_continuous(breaks = seq(0, 100, 20),
                       sec.axis = sec_axis(~ . , name = "Embryo size",breaks = NULL, labels = NULL)) +
    coord_cartesian(ylim = c(0, 100)) + 
    labs(y = "Percentage of embryos\nmatching biopsy class",
         x = "Biopsy aneuploidy") +
    facet_grid(Embryo_size~Dispersal) +
    theme_classic() +
    theme(legend.position = "none")

  save.double.width(col.plot.pgdis, 
    filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_columns_two_biopsy_pgdis_", b), 
    height = 150)
  
  col.data.merge = in.data %>% 
    dplyr::filter(Biopsy_size == b) %>%
    dplyr::mutate(EmbryoMergeClass = to.merged.class(Aneuploidy),
                  BiopsyMergeClass = to.merged.class(n_aneuploid/Biopsy_size)) %>%
    dplyr::group_by(Dispersal, n_aneuploid, Embryo_size, EmbryoMergeClass, BiopsyMergeClass) %>%
    dplyr::summarise(Count = sum(Total_biopsies)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Embryo_size, n_aneuploid, Dispersal) %>%
    dplyr::mutate(Total = sum(Count),
                  PctTotal = Count / Total * 100,
                  IsPGDISCorrect = EmbryoMergeClass == BiopsyMergeClass) %>%
    dplyr::filter(IsPGDISCorrect)
  
  col.plot.merge = ggplot(col.data.merge, aes(x = n_aneuploid/b*100, y= PctTotal,  fill=IsPGDISCorrect))+
    geom_hline(yintercept = 50) +
    geom_col(position="stack") +
    scale_fill_manual(values = c("dark green")) +
    scale_x_continuous(breaks = seq(0, 100, 20), 
                       sec.axis = sec_axis(~ . , name = "Dispersal of aneuploid cells",breaks = NULL, labels = NULL)) +
    scale_y_continuous(breaks = seq(0, 100, 20),
                       sec.axis = sec_axis(~ . , name = "Embryo size",breaks = NULL, labels = NULL)) +
    coord_cartesian(ylim = c(0, 100)) + 
    labs(y = "Percentage of embryos\nmatching biopsy class",
         x = "Biopsy aneuploidy") +
    facet_grid(Embryo_size~Dispersal) +
    theme_classic() +
    theme(legend.position = "none")

  save.double.width(col.plot.merge, 
    filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_columns_two_biopsy_merge_", b), 
    height = 150)
  
}
