# Make figures showing the impact of embryo and biopsy parameters.
library(tessera)
library(tidyverse)
library(patchwork)
library(svglite)
library(data.table)
library(parallel)

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




################################################################################
# Figure 4: all biopsy heatmaps. Generated while looping over
# biopsy sizes

################################################################################
# Figure 5: biopsy origin heatmap. Generated by analyseCombos while looping over
# biopsy sizes

################################################################################
# Figure 6: predictive heatmap and columns

# Calculate averages
aneuploidy.dispersal.biopsy.values = agg.values %>% 
  group_by(Aneuploidy, Dispersal, Biopsy_size) %>%
  summarise(pct_pgdis_match = mean(f_pgdis_match)*100,
            sd_pgdis_match = sd(f_pgdis_match)*100,
            pct_merge_match = mean(f_merge_match)*100,
            sd_merge_match = sd(f_merge_match)*100)


# Make a heatmap for 5 cell biopsies
heatmap_pgdis_match_5 = ggplot(aneuploidy.dispersal.biopsy.values[aneuploidy.dispersal.biopsy.values$Biopsy_size==5,], 
       aes(x = Aneuploidy*100, y = Dispersal, fill=pct_pgdis_match)) +
  geom_raster()+
  geom_vline(xintercept = 19.5, col="white")+
  geom_vline(xintercept = 39.5, col="white")+
  geom_vline(xintercept = 80.5, col="white")+
  annotate(geom = "text", label = "Euploid", x=10, y=1.06)+
  annotate(geom = "text", label = "Low", x=30, y=1.06)+
  annotate(geom = "text", label = "High", x=60, y=1.06)+
  annotate(geom = "text", label = "Aneuploid", x=91, y=1.06)+
  labs(x = "Aneuploidy",
       y = "Dispersal",
       fill = "Percent of\nbiopsies in\ncorrect\nclass") + 
  coord_cartesian(ylim = c(0, 1.05))+
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_fill_viridis_c()+
  theme_classic()
save.single.width(heatmap_pgdis_match_5, 
                  filename = paste0(FIGURE.OUTPUT.DIR, "heatmap_pgdis_match_5"), 
                  height = 60)

heatmap_merge_match_5 = ggplot(aneuploidy.dispersal.biopsy.values[aneuploidy.dispersal.biopsy.values$Biopsy_size==5,], 
       aes(x = Aneuploidy*100, y = Dispersal, fill=pct_merge_match)) +
  geom_raster()+
  geom_vline(xintercept = 19.5, col="white")+
  geom_vline(xintercept = 80.5, col="white")+
  annotate(geom = "text", label = "Euploid", x=10, y=1.06)+
  annotate(geom = "text", label = "Mosaic", x=50, y=1.06)+
  annotate(geom = "text", label = "Aneuploid", x=91, y=1.06)+
  labs(x = "Aneuploidy",
       y = "Dispersal",
       fill = "Percent of\nbiopsies in\ncorrect\nclass") + 
  coord_cartesian(ylim = c(0, 1.05))+
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_fill_viridis_c()+
  theme_classic()
save.double.width(heatmap_merge_match_5, 
                  filename = paste0(FIGURE.OUTPUT.DIR, "heatmap_merge_match_5.png"), 
                  height = 120)

# Using the PGDIS classifictions
aneuploidy.dispersal.biopsy.heatmap = ggplot(aneuploidy.dispersal.biopsy.values, aes(x = Aneuploidy*100, y = Dispersal, fill=pct_pgdis_match)) +
  geom_raster()+
  geom_vline(xintercept = 19.5, col="white")+
  geom_vline(xintercept = 39.5, col="white")+
  geom_vline(xintercept = 80.5, col="white")+
  labs(x = "Aneuploidy",
       y = "Dispersal",
       fill = "Percent of\nbiopsies in\ncorrect\nclass") + 
  coord_cartesian(ylim = c(0, 1))+
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_fill_viridis_c()+
  facet_wrap(~Biopsy_size, ncol = 4)+
  theme_classic()

save.double.width(aneuploidy.dispersal.biopsy.heatmap, 
                  filename = paste0(FIGURE.OUTPUT.DIR, "heatmap_pgdis_match"), 
                  height = 120)

# As above, but using the merged classifictions instead of PGDIS
aneuploidy.dispersal.biopsy.merge.heatmap = ggplot(aneuploidy.dispersal.biopsy.values, 
                                                   aes(x = Aneuploidy*100, y = Dispersal, fill=pct_merge_match)) +
  geom_raster()+
  geom_vline(xintercept = 19.5, col="white")+
  geom_vline(xintercept = 80.5, col="white")+
  labs(x = "Aneuploidy",
       y = "Dispersal",
       fill = "Percent of\nbiopsies in\ncorrect\nclass") + 
  coord_cartesian(ylim = c(0, 1))+
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_fill_viridis_c()+
  facet_wrap(~Biopsy_size, ncol = 4)+
  theme_classic()

save.double.width(aneuploidy.dispersal.biopsy.merge.heatmap, 
                  filename = paste0(FIGURE.OUTPUT.DIR, "heatmap_merge_match.png"), 
                  height = 120)