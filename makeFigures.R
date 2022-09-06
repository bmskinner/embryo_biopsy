# Make other figures for the paper not created in other analysis files
library(tessera)
library(tidyverse)
library(patchwork)
library(svglite)
library(data.table)
library(parallel)

# source("analyseCombos.R")

save.single.width = function(plot, filename, height){
  ggsave(plot, filename=paste0(filename, ".png"), dpi=300, units = "mm", width = 85, height = height)
  ggsave(plot, filename=paste0(filename, ".svg"), dpi=300, units = "mm", width = 85, height = height)
}
save.double.width = function(plot, filename, height){
  ggsave(plot, filename=paste0(filename, ".png"), dpi=300, units = "mm", width = 170, height = height)
  ggsave(plot, filename=paste0(filename, ".svg"), dpi=300, units = "mm", width = 170, height = height)
}

# Make a basic heatmap with no rectangle annotations
plot.heatmap = function(data, zero.data){
  ggplot(data, aes(x = f_aneuploid, y = Aneuploidy, fill=PctTotal))+
    geom_raster(data = zero.data)+
    geom_raster() +
    scale_fill_viridis_c() +
    labs(x = "Biopsy aneuploidy",
         y = "Embryo aneuploidy",
         fill = "Percentage\nof biopsies") +
    scale_y_continuous(breaks = seq(0, 1, 0.1)) +
    scale_x_continuous(breaks = seq(0, 100, 20)) +
    theme_classic() + 
    facet_wrap(~Dispersal, ncol = 3)+
    theme(legend.position = "top",
          legend.title = element_text(size=9),
          legend.text = element_text(size=9),
          legend.box.spacing = unit(1.5, "mm"),
          legend.key.height = unit(3, "mm"),
          legend.title.align = 1,
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
}

# Make a heatmap with PGDIS class annotations
plot.pgdis.heatmap = function(data, zero.data){
  plot.heatmap(data, zero.data)+
    geom_rect(xmin=-10, xmax=10, ymin=-0.025, ymax=0.195, 
              fill=NA, col="white", size=1)+
    geom_rect(xmin=10, xmax=30, ymin=0.195, ymax=0.395, 
              fill=NA, col="white", size=1)+
    geom_rect(xmin=30, xmax=90, ymin=0.395, ymax=0.805, 
              fill=NA, col="white", size=1)+
    geom_rect(xmin=90, xmax=110, ymin=0.805, ymax=1.025, 
              fill=NA, col="white", size=1)
}

# Make a heatmap with merged class annotations 20% and 80%
plot.merge.heatmap = function(data, zero.data){
  plot.heatmap(data, zero.data)+
    geom_rect(xmin=-10, xmax=10, ymin=-0.025, ymax=0.195,
              fill=NA, col="white", size=1)+
    geom_rect(xmin=10, xmax=90, ymin=0.195, ymax=0.825,
              fill=NA, col="white", size=1)+
    geom_rect(xmin=90, xmax=110, ymin=0.825, ymax=1.025,
              fill=NA, col="white", size=1)
}

# Make a heatmap with merged class annotations 10% and 90%
plot.merge2.heatmap = function(data, zero.data){
  plot.heatmap(data, zero.data)+
    geom_rect(xmin=-10, xmax=10, ymin=-0.025, ymax=0.095,
              fill=NA, col="white", size=1)+
    geom_rect(xmin=10, xmax=90, ymin=0.095, ymax=0.905,
              fill=NA, col="white", size=1)+
    geom_rect(xmin=90, xmax=110, ymin=0.905, ymax=1.025,
              fill=NA, col="white", size=1)
}


plot.columns = function(data){
  ggplot(data, aes(x = n_aneuploid/b*100, y= PctTotal,  fill=IsCorrect))+
    geom_hline(yintercept = 50) +
    geom_col(position="stack") +
    scale_fill_manual(values = c("dark green")) +
    scale_x_continuous(breaks = seq(0, 100, 20)) +
    scale_y_continuous(breaks = seq(0, 100, 20)) +
    coord_cartesian(ylim = c(0, 100)) + 
    labs(y = "Percentage of embryos\nmatching biopsy class",
         x = "Biopsy aneuploidy") +
    facet_wrap(~Dispersal, ncol = 3) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}



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


biopsy.values = do.call(rbind, mclapply(list.files(path = "data/aggregates", 
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
# # Read the saved raw values for selected dispersals
raw.values = do.call(rbind, mclapply(list.files(path = "data/raw",
                                                pattern = "raw_values_a.*_d(0|0.5|1).csv", full.names = T),
                                     fread,
                                     header = T,
                                     mc.cores = 5))

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
  
  ggsave(hmap.pgdis.plot, filename = paste0("figure/predictive_heatmap_pgdis_",b,".svg"), dpi=300, units = "mm", width = 85, height = 85)
  ggsave(hmap.pgdis.plot, filename = paste0("figure/predictive_heatmap_pgdis_",b,".png"), dpi=300, units = "mm", width = 85, height = 85)
  
  ggsave(hmap.merge.plot, filename = paste0("figure/predictive_heatmap_merge_",b,".svg"), dpi=300, units = "mm", width = 85, height = 85)
  ggsave(hmap.merge.plot, filename = paste0("figure/predictive_heatmap_merge_",b,".png"), dpi=300, units = "mm", width = 85, height = 85)
  
  ggsave(hmap.merge2.plot, filename = paste0("figure/predictive_heatmap_merge2_",b,".svg"), dpi=300, units = "mm", width = 85, height = 85)
  ggsave(hmap.merge2.plot, filename = paste0("figure/predictive_heatmap_merge2_",b,".png"), dpi=300, units = "mm", width = 85, height = 85)
  
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
  ggsave(col.pgdis.plot, filename = paste0("figure/predictive_columns_pgdis_",b,".svg"), dpi=300, units = "mm", width = 85, height = 85)
  ggsave(col.pgdis.plot, filename = paste0("figure/predictive_columns_pgdis_",b,".png"), dpi=300, units = "mm", width = 85, height = 85)
  
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
  ggsave(col.merge.plot, filename = paste0("figure/predictive_columns_merge_",b,".svg"), dpi=300, units = "mm", width = 85, height = 85)
  ggsave(col.merge.plot, filename = paste0("figure/predictive_columns_merge_",b,".png"), dpi=300, units = "mm", width = 85, height = 85)
  
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
  ggsave(col.merge2.plot, filename = paste0("figure/predictive_columns_merge2_",b,".svg"), dpi=300, units = "mm", width = 85, height = 85)
  ggsave(col.merge2.plot, filename = paste0("figure/predictive_columns_merge2_",b,".png"), dpi=300, units = "mm", width = 85, height = 85)
  
  combo.pgdis.plot = hmap.pgdis.plot + col.pgdis.plot + plot_annotation(tag_levels = c("A"))
  ggsave(combo.pgdis.plot, filename = paste0("figure/predictive_heatmap_pgdis_combined_",b,".svg"), dpi=300, units = "mm", 
         width = 170, height = 85)
  ggsave(combo.pgdis.plot, filename = paste0("figure/predictive_heatmap_pgdis_combined_",b,".png"), dpi=300, units = "mm", 
         width = 170, height = 85)
  
  combo.merge.plot = hmap.merge.plot + col.merge.plot + plot_annotation(tag_levels = c("A"))
  ggsave(combo.merge.plot, filename = paste0("figure/predictive_heatmap_merge_combined_",b,".svg"), dpi=300, units = "mm", 
         width = 170, height = 85)
  ggsave(combo.merge.plot, filename = paste0("figure/predictive_heatmap_merge_combined_",b,".png"), dpi=300, units = "mm", 
         width = 170, height = 85)
  
  combo.merge2.plot = hmap.merge2.plot + col.merge2.plot + plot_annotation(tag_levels = c("A"))
  ggsave(combo.merge2.plot, filename = paste0("figure/predictive_heatmap_merge2_combined_",b,".svg"), dpi=300, units = "mm", 
         width = 170, height = 85)
  ggsave(combo.merge2.plot, filename = paste0("figure/predictive_heatmap_merge2_combined_",b,".png"), dpi=300, units = "mm", 
         width = 170, height = 85)
  
}

# Calculate step values for aneuploidy
# All values for filling out the steps
all.vals = expand.grid("Aneuploidy" = ANEUPLOIDY.RANGE,
                       "Dispersal" = 1,
                       "Biopsy_size" = 5,
                       "diff_to_embryo" = round(seq(0, 1, 0.01), digits = 2),
                       "Count" = 0)

step.aneuploidy = calc.step.accuracy(raw.values[raw.values$Dispersal==1 & Biopsy_size == 5,], all.vals)

# Calculate step values for dispersal
# All values for filling out the steps
all.vals = expand.grid("Aneuploidy" = 0.2,
                       "Dispersal" = DISPERSAL.RANGE,
                       "Biopsy_size" = 5,
                       "diff_to_embryo" = round(seq(0, 1, 0.01), digits = 2),
                       "Count" = 0)

# Read in the files needed for dispersal step calculation
step.values = do.call(rbind, mclapply(list.files(path = "data/raw",
                                                 pattern = "raw_values_a0.2_d.*.csv", full.names = T),
                                      fread,
                                      header = T,
                                      mc.cores = 5))
step.dispersal  = calc.step.accuracy(step.values[step.values$Aneuploidy == 0.2 & Biopsy_size == 5,], all.vals)



################################################################################
# Figure 1: embryo models and their biopsy distribution
# Dispersed embryo
dispersed = tessera::Embryo(n.cells = 200, 
                    n.chrs = 1,
                    prop.aneuploid = 0.2, 
                    dispersal = 1, 
                    concordance = 1, 
                    rng.seed = 42)
plot(dispersed)
# Clustered embryo
clustered = tessera::Embryo(n.cells = 200, 
                            n.chrs = 1,
                            prop.aneuploid = 0.2, 
                            dispersal = 0, 
                            concordance = 1, 
                            rng.seed = 42)
plot(clustered)

dispersed.biopsies = tessera::takeAllBiopsies(dispersed, biopsy.size = 5, chromosome = 1)
clustered.biopsies = tessera::takeAllBiopsies(clustered, biopsy.size = 5, chromosome = 1)

p = plot.biopsy.aneuploidy(dispersed.biopsies, 5, 0.2) + plot.biopsy.aneuploidy(clustered.biopsies, 5, 0.2)
save.single.width(p, filename="figure/Figure_1CD", height = 40)

################################################################################
# Figure 2: Step plots for aneuploidy and dispersal
# Step values calculated in analyseCombos.R
# Filter down so the plot is not overcrowded
step.aneuploidy.filt = step.aneuploidy[step.aneuploidy$Aneuploidy %in% seq(0, 0.5, 0.05),]

aneuploidy.step.plot = ggplot(step.aneuploidy.filt, aes(x = diff_to_embryo*100, y = CumPct, 
                                                        col=Aneuploidy*100, group=Aneuploidy)) +
  geom_step(size=2)+
  labs(y = "Percent of biopsies" ,
       x = "Percent difference\nbetween biopsy and embryo",
       col = "Aneuploidy") +
  scale_colour_viridis_c() +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  theme_classic()+
  theme(legend.position = c(0.85, 0.4),
        legend.title = element_text(size=9),
        axis.title = element_text(size=9))

# Supplementary figure: aneuploidy step plot 0.5-1
step.aneuploidy.supp = step.aneuploidy[step.aneuploidy$Aneuploidy %in% seq(0.5, 1.0, 0.05),]
aneuploidy.step.plot.supp = ggplot(step.aneuploidy.supp, aes(x = diff_to_embryo*100, y = CumPct, 
                                 col=Aneuploidy*100, group=Aneuploidy)) +
  geom_step(size=2)+
  labs(y = "Percent of biopsies" ,
       x = "Percent difference\nbetween biopsy and embryo",
       col = "Aneuploidy") +
  scale_colour_viridis_c() +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  theme_classic()+
  theme(legend.position = c(0.85, 0.4),
        legend.title = element_text(size=9),
        axis.title = element_text(size=9))
save.single.width(aneuploidy.step.plot.supp, filename="figure/Figure_S1_step", height = 85)


# Dispersal step plot

dispersal.step.plot = ggplot(step.dispersal, aes(x = diff_to_embryo*100, y = CumPct, 
                                                        col=Dispersal, group=Dispersal)) +
  geom_step(size=2)+
  labs(y = "Percent of biopsies" ,
       x = "Percent difference\nbetween biopsy and embryo") +
  scale_colour_viridis_c() +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  theme_classic()+
  theme(legend.position = c(0.85, 0.4),
        legend.title = element_text(size=9),
        axis.title = element_text(size=9))

fig2 = aneuploidy.step.plot + dispersal.step.plot + plot_annotation(tag_levels = c("A"))
save.double.width(fig2, filename="figure/Figure_2_step", height = 85)

################################################################################
# Figure 3: 5-cell biopsy heatmap. Generated by analyseCombos while looping over
# biopsy sizes

# Read in the aggregate values
agg.values = do.call(rbind, mclapply(list.files(path = "data/aggregates", 
                                                pattern = "merged", full.names = T), 
                                     fread, 
                                     header = T, mc.cores = N.CORES))

# Filter to 5-cell biopsy and calculate matching percentage
aneuploidy.dispersal.values = agg.values %>% 
  dplyr::filter(Biopsy_size == 5) %>%
  dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size) %>%
  dplyr::summarise(pct_pgdis_match = mean(f_pgdis_match)*100,
                   sd_pgdis_match = sd(f_pgdis_match)*100)

aneuploidy.dispersal.heatmap = ggplot(aneuploidy.dispersal.values, 
                                      aes(x = Aneuploidy*100, y = Dispersal, fill=pct_pgdis_match)) +
  geom_raster()+
  geom_vline(xintercept = 19.5, col="white")+
  geom_vline(xintercept = 39.5, col="white")+
  geom_vline(xintercept = 80.5, col="white")+
  annotate(geom = "text", label = "Euploid", x=10, y=1.06, size=3)+
  annotate(geom = "text", label = "Low", x=30, y=1.06, size=3)+
  annotate(geom = "text", label = "High", x=60, y=1.06, size=3)+
  annotate(geom = "text", label = "Aneuploid", x=91, y=1.06, size=3)+
  labs(x = "Embryo aneuploidy",
       y = "Embryo dispersal",
       fill = "Percent of biopsies\nin correct class") + 
  coord_cartesian(ylim = c(0, 1.05))+
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_fill_viridis_c(limits = c(0, 100))+
  theme_classic()+
  theme(legend.position = "top",
        legend.title = element_text(size=9),
        legend.text = element_text(size=9),
        legend.box.spacing = unit(1.5, "mm"),
        legend.key.height = unit(3, "mm"),
        legend.title.align = 1,
        axis.title = element_text(size=9))
save.single.width(aneuploidy.dispersal.heatmap, filename="figure/Figure_3_heatmap", height = 85)

################################################################################
# Figure 4: all biopsy heatmaps. Generated by analyseCombos while looping over
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
ggplot(aneuploidy.dispersal.biopsy.values[aneuploidy.dispersal.biopsy.values$Biopsy_size==5,], 
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
ggsave(filename = "figure/heatmap_pgdis_match_5.png", dpi=300, units = "mm", width = 85, height = 60)

ggplot(aneuploidy.dispersal.biopsy.values[aneuploidy.dispersal.biopsy.values$Biopsy_size==5,], 
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
ggsave(filename = "figure/heatmap_merge_match_5.png", dpi=300, units = "mm", width = 170, height = 120)


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
ggsave(aneuploidy.dispersal.biopsy.heatmap, filename = "figure/heatmap_pgdis_match.png", dpi=300, units = "mm", width = 170, height = 120)

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
ggsave(aneuploidy.dispersal.biopsy.merge.heatmap, filename = "figure/heatmap_merge_match.png", dpi=300, units = "mm", width = 170, height = 120)
# saveRDS(aneuploidy.dispersal.biopsy.heatmap, "data/aneuploidy.dispersal.biopsy.heatmap.Rds")



