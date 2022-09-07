# Make figures showing the impact of embryo and biopsy parameters.
library(tessera)
library(tidyverse)
library(patchwork)
library(svglite)
library(data.table)
library(parallel)
library(fs)

source("parameters.R")


# Plot the percentage of cells in a biopsy versus the percentage of total biopsies
plot.biopsy.aneuploidy = function(n.biopsy.aneuploid, biopsy.size, pct.embryo.aneuploid){
  pct.biopsy.aneuploid = n.biopsy.aneuploid / biopsy.size * 100
  biopsy.data = as.data.frame(pct.biopsy.aneuploid) %>%
    dplyr::group_by(pct.biopsy.aneuploid) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::mutate(PctBiopsies = n/sum(n)*100)
  
  # Number of biopsies
  ggplot(biopsy.data, aes(x=pct.biopsy.aneuploid, y = PctBiopsies)) +
    geom_col(width = 5) +
    coord_cartesian(xlim = c(0, 100), 
                    ylim = c(0, 100)) +
    scale_x_continuous(breaks = seq(0, 100, 20)) +
    scale_y_continuous(breaks = seq(0, 100, 20)) + 
    labs(x = "Percent aneuploid cells in biopsy",
         y = "Percent of biopsies") +
    theme_classic()
}

if(!fs::dir_exists(FIGURE.OUTPUT.DIR)){
  fs::dir_create(FIGURE.OUTPUT.DIR, recursive = T)
}

biopsy.values = do.call(rbind, mclapply(list.files(path = AGGREGATE.DATA.PATH, 
                                                   pattern = "biopsy", full.names = T), 
                                        fread, 
                                        header = T, mc.cores = N.CORES))

# Calculate the proportion of biopsies originating from each embryo combination
for(b in BIOPSY.SIZES){
  
  biopsy.aggregate = biopsy.values %>% 
    dplyr::filter(Biopsy_size == b) %>%
    dplyr::select(-embryo_merge_class, -biopsy_merge_class, -n_biopsies_merge, 
                  -embryo_pgdis_class, -biopsy_pgdis_class) %>%
    dplyr::group_by(Biopsy_size, n_aneuploid) %>%
    dplyr::mutate(total_biopsies = sum(n_biopsies_pgdis),
                  pct_biopsies = n_biopsies_pgdis/total_biopsies*100)
  
  zero.data = expand.grid(Aneuploidy = ANEUPLOIDY.RANGE,
                          Dispersal = DISPERSAL.RANGE,
                          pct_biopsies = 0)
  
  ggplot(biopsy.aggregate, aes(x=Aneuploidy, y = Dispersal, fill = pct_biopsies))+
    geom_raster(data = zero.data)+
    geom_raster()+
    labs(fill = "Percent of\nbiopsies")+
    scale_fill_viridis_c()+
    facet_wrap(~n_aneuploid)+
    theme_classic()
  ggsave(filename = paste0("figure/heatmap_biopsy_",b,".png"), dpi=300, units = "mm", width = 170, height = 120)
}

# Make the predictive heatmaps
# # Read the saved raw values for selected dispersals at 200 cell embryo
raw.values = do.call(rbind, mclapply(list.files(path = RAW.DATA.PATH,
                                                pattern = "raw_values_e200_a.*_d(0|0.5|1).csv", full.names = T),
                                     fread,
                                     header = T,
                                     mc.cores = N.CORES))

# Create predictive heatmaps and column charts for PGDIS and merged classes
for(b in BIOPSY.SIZES){
  # Calculate the proportion of embryos matching each biopsy
  agg.values.heatmap.data = raw.values %>% 
    dplyr::filter(Biopsy_size == b) %>%
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
                          f_aneuploid = seq(0, 100, 100/b),
                          PctTotal = 0)
  
  # Plot with PGDIS classes
  hmap.pgdis.plot  = plot.pgdis.heatmap(agg.values.heatmap.data, zero.data)
  hmap.merge.plot  = plot.merge.heatmap(agg.values.heatmap.data, zero.data)
  hmap.merge2.plot = plot.merge2.heatmap(agg.values.heatmap.data, zero.data)
  
  ggsave(hmap.pgdis.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_heatmap_pgdis_",b,".svg"), dpi=300, units = "mm", width = 85, height = 85)
  ggsave(hmap.pgdis.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_heatmap_pgdis_",b,".png"), dpi=300, units = "mm", width = 85, height = 85)
  
  ggsave(hmap.merge.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_heatmap_merge_",b,".svg"), dpi=300, units = "mm", width = 85, height = 85)
  ggsave(hmap.merge.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_heatmap_merge_",b,".png"), dpi=300, units = "mm", width = 85, height = 85)
  
  ggsave(hmap.merge2.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_heatmap_merge2_",b,".svg"), dpi=300, units = "mm", width = 85, height = 85)
  ggsave(hmap.merge2.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_heatmap_merge2_",b,".png"), dpi=300, units = "mm", width = 85, height = 85)
  
  # Calculate the proportion of embryos matching biopsy class
  agg.values.col.data = raw.values %>% 
    dplyr::filter(Biopsy_size == b) %>%
    rowwise %>%
    mutate(n_aneuploid = list(as.double(c_across(starts_with("V"))))) %>%
    select(-starts_with("V")) %>%
    tidyr::unnest(n_aneuploid) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(EmbryoClass = to.pgdis.class(Aneuploidy),
                  BiopsyClass = to.pgdis.class(n_aneuploid/Biopsy_size)) %>%
    dplyr::group_by(Dispersal, n_aneuploid, EmbryoClass, BiopsyClass) %>%
    dplyr::summarise(Count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(n_aneuploid, Dispersal) %>%
    dplyr::mutate(Total = sum(Count),
                  PctTotal = Count / Total * 100,
                  IsCorrect = EmbryoClass == BiopsyClass) %>%
    dplyr::filter(IsCorrect)
  
  col.pgdis.plot = plot.columns(agg.values.col.data)
  ggsave(col.pgdis.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_columns_pgdis_",b,".svg"), dpi=300, units = "mm", width = 85, height = 85)
  ggsave(col.pgdis.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_columns_pgdis_",b,".png"), dpi=300, units = "mm", width = 85, height = 85)
  
  agg.values.merge.data = raw.values %>% 
    dplyr::filter(Biopsy_size == b) %>%
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
  
  col.merge.plot = plot.columns(agg.values.merge.data)
  ggsave(col.merge.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_columns_merge_",b,".svg"), dpi=300, units = "mm", width = 85, height = 85)
  ggsave(col.merge.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_columns_merge_",b,".png"), dpi=300, units = "mm", width = 85, height = 85)
  
  agg.values.merge2.data = raw.values %>% 
    dplyr::filter(Biopsy_size == b) %>%
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
  
  col.merge2.plot = plot.columns(agg.values.merge2.data)
  ggsave(col.merge2.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_columns_merge2_",b,".svg"), dpi=300, units = "mm", width = 85, height = 85)
  ggsave(col.merge2.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_columns_merge2_",b,".png"), dpi=300, units = "mm", width = 85, height = 85)
  
  combo.pgdis.plot = hmap.pgdis.plot + col.pgdis.plot + plot_annotation(tag_levels = c("A"))
  ggsave(combo.pgdis.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_heatmap_pgdis_combined_",b,".svg"), dpi=300, units = "mm", 
         width = 170, height = 85)
  ggsave(combo.pgdis.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_heatmap_pgdis_combined_",b,".png"), dpi=300, units = "mm", 
         width = 170, height = 85)
  
  combo.merge.plot = hmap.merge.plot + col.merge.plot + plot_annotation(tag_levels = c("A"))
  ggsave(combo.merge.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_heatmap_merge_combined_",b,".svg"), dpi=300, units = "mm", 
         width = 170, height = 85)
  ggsave(combo.merge.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_heatmap_merge_combined_",b,".png"), dpi=300, units = "mm", 
         width = 170, height = 85)
  
  combo.merge2.plot = hmap.merge2.plot + col.merge2.plot + plot_annotation(tag_levels = c("A"))
  ggsave(combo.merge2.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_heatmap_merge2_combined_",b,".svg"), dpi=300, units = "mm", 
         width = 170, height = 85)
  ggsave(combo.merge2.plot, filename = paste0(FIGURE.OUTPUT.DIR, "/predictive_heatmap_merge2_combined_",b,".png"), dpi=300, units = "mm", 
         width = 170, height = 85)
  
}
