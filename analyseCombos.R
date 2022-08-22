# Test analysing pre-computed values
library(data.table)
library(tidyverse)
library(parallel)
library(svglite)
library(patchwork)

ANEUPLOIDY.RANGE = seq(0, 1, 0.01)
DISPERSAL.RANGE = seq(0, 1, 0.01)
N.REPLICATES = 100 
BIOPSY.SIZES = c(3:10, 15, 20, 25, 30)
N.CORES = ifelse(Sys.info()["sysname"]=="Windows", 1, 20) 

to.pgdis.class = function(f.aneuploidy){
  case_when(f.aneuploidy < 0.20 ~ "Euploid",
            f.aneuploidy < 0.40 ~ "Low level",
            f.aneuploidy <= 0.80 ~ "High level",
            f.aneuploidy <= 1 ~ "Aneuploid")
}

to.merged.class = function(f.aneuploidy){
  case_when(f.aneuploidy < 0.20 ~ "Euploid",
            f.aneuploidy <= 0.80 ~ "Mosaic",
            f.aneuploidy <= 1 ~ "Aneuploid")
}

to.merged.class.2 = function(f.aneuploidy){
  case_when(f.aneuploidy < 0.10 ~ "Euploid",
            f.aneuploidy <= 0.90 ~ "Mosaic",
            f.aneuploidy <= 1 ~ "Aneuploid")
}

# data.file = "data/all.combos.csv"

# Given raw input data, calculate the difference between the biopsy and 
# embryo aneuploidies and summarise into counts
calc.step.accuracy = function(data){

  data %>%
    rowwise %>%
    mutate(f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size))) %>%
    select(-starts_with("V")) %>%
    ungroup %>%
    unnest(f_aneuploid) %>%
    dplyr::mutate(diff_to_embryo = round(abs(f_aneuploid-Aneuploidy), digits=2)) %>%
    dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size, diff_to_embryo) %>%
    dplyr::summarise(Count = dplyr::n()) %>%
    dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size) %>%
    dplyr::mutate(TotalBiopsies = sum(Count),
                  PctBiopsies   = Count/TotalBiopsies*100) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Dispersal, Biopsy_size, Aneuploidy, diff_to_embryo) %>%
    dplyr::select(-TotalBiopsies, -Count) %>%
    dplyr::distinct()
  
}

# Given raw input data, calculate the accuracy of the biopsies
# with respect to PGDIS classification boundaries
calc.pgdis.accuracy = function(data){
  

  
  data %>%
    rowwise %>%
    mutate(n_aneuploid = list(c_across(starts_with("V"))),
           f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size))) %>%
    select(-starts_with("V")) %>%
    mutate(pgdis_class = list(to.pgdis.class(f_aneuploid)),
           merge_class = list(to.merged.class(f_aneuploid)),
           actual_pgdis_class = to.pgdis.class(Aneuploidy),
           actual_merge_class = to.merged.class(Aneuploidy),
           n_pgdis_match = sum(pgdis_class == actual_pgdis_class),
           f_pgdis_match = n_pgdis_match / 200,
           n_merge_match = sum(merge_class == actual_merge_class),
           f_merge_match = n_merge_match / 200) %>%
    select(-n_aneuploid, -f_aneuploid, -pgdis_class, -merge_class) %>%
    ungroup
}

# Calculate the number of biopsies in each class for each aneuploidy, dispersal
# and biopsy size combination
calc.biopsy.accuracy = function(data){
  
  data %>%
    rowwise %>%
    dplyr::mutate(n_aneuploid = list(c_across(starts_with("V")))) %>%
    dplyr::select(-starts_with("V")) %>%
    tidyr::unnest(n_aneuploid) %>%
    dplyr::mutate(embryo_pgdis_class = to.pgdis.class(Aneuploidy),
                  embryo_merge_class = to.merged.class(Aneuploidy),
                  biopsy_pgdis_class = to.pgdis.class(n_aneuploid/Biopsy_size),
                  biopsy_merge_class = to.merged.class(n_aneuploid/Biopsy_size)) %>%
    dplyr::ungroup() %>% # cancels rowwise
    dplyr::group_by(Aneuploidy, Dispersal, n_aneuploid, embryo_pgdis_class, biopsy_pgdis_class, Biopsy_size) %>%
    dplyr::mutate(n_biopsies_pgdis = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Aneuploidy, Dispersal, n_aneuploid, embryo_merge_class, biopsy_merge_class, Biopsy_size) %>%
    dplyr::mutate(n_biopsies_merge = dplyr::n()) %>%
    dplyr::select(-Seed) %>%
    dplyr::distinct()
}

# Calculate aggregate data from raw values for making heatmaps
make.aggregate.values = function(a, d){
  aneu.part.file = paste0("data/aggregates/merged_a", a, "_d", d, ".csv")
  if(!file.exists(aneu.part.file)){
    in.file = paste0("data/raw/raw_values_a", a, "_d", d, ".csv")
    if(file.exists(in.file)){
      in.data = fread(file=in.file, header = T)
      filt.tf = calc.pgdis.accuracy(in.data)
      write.csv(filt.tf, file = aneu.part.file, quote = F, row.names = F)
      rm(filt.tf)
      gc()
    }
  }
}

# Make output files aggregating data for biopsy accuracy calculations
make.biopsy.values = function(a){
  # out.file = paste0("data/aggregates/biopsy_a", a, "_d", d, ".csv")
  out.file = paste0("data/aggregates/biopsy_a", a, ".csv")
  if(!file.exists(out.file)){
    in.files = list.files(path = "data/raw",
                          pattern = paste0("raw_values_a", a, "_"), full.names = T)
    
    in.data = do.call(rbind, mclapply(in.files,
                                      fread,
                                      header = T,
                                      mc.cores = 22))
    filt.tf = calc.biopsy.accuracy(in.data)
    write.csv(filt.tf, file = out.file, quote = F, row.names = F)
    rm(filt.tf)
    invisible(gc(verbose = F))
  }
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
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
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
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

combinations = expand.grid(a = ANEUPLOIDY.RANGE, d = DISPERSAL.RANGE)

# Functions write output files, no need to store in object
mcmapply(make.aggregate.values, a = combinations$a, d = combinations$d, mc.cores = N.CORES)
mcmapply(make.biopsy.values, a = ANEUPLOIDY.RANGE, mc.cores = 3)

# Now read in the aggregate values to generate the summary figures
agg.values = do.call(rbind, mclapply(list.files(path = "data/aggregates", 
                                              pattern = "merged", full.names = T), 
                                   fread, 
                                   header = T, mc.cores = N.CORES))


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
  ggsave(col.pgdis.plot, filename = paste0("figure/predictive_columns_merge_",b,".svg"), dpi=300, units = "mm", width = 85, height = 85)
  ggsave(col.pgdis.plot, filename = paste0("figure/predictive_columns_merge_",b,".png"), dpi=300, units = "mm", width = 85, height = 85)
  
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
  ggsave(col.pgdis.plot, filename = paste0("figure/predictive_columns_merge2_",b,".svg"), dpi=300, units = "mm", width = 85, height = 85)
  ggsave(col.pgdis.plot, filename = paste0("figure/predictive_columns_merge2_",b,".png"), dpi=300, units = "mm", width = 85, height = 85)
  
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
# These plots are hard to interpret - keep for supplement
step.aneuploidy = calc.step.accuracy(raw.values[raw.values$Dispersal==1 & Biopsy_size == 5,])

step.aneuploidy = step.aneuploidy %>% mutate(Bin = cut(diff_to_embryo, breaks=seq(-.01, 1.01, 0.1)))

# Make a table of all possible values for cumulative summing
# all.vals = expand.grid("Aneuploidy" = ANEUPLOIDY.RANGE,
#                        "Dispersal" = unique(grouped_data$Dispersal),
#                        "Biopsy_size" = unique(grouped_data$Biopsy_size),
#                        "diff_to_embryo" = round(seq(0, 1, 0.01), digits = 2),
#                        "PctBiopsies" = 0)

ggplot(step.aneuploidy, aes(x = Aneuploidy, y = PctBiopsies, fill = Bin)) +
  geom_col(position = "stack", width=0.01)+
  labs(y = "Percent of biopsies" ,
       x = "Aneuploidy of embryo",
       fill = "Percent difference\nbetween biopsy\nand embryo") +
  scale_fill_viridis_d(direction = -1) +
  # coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  # scale_x_continuous(breaks = seq(0, 100, 10)) +
  # scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_classic()

# aneuploidy.values = agg.values %>%
#   filter(Dispersal %in% c(1),
#          Biopsy_size==5) %>%
#   group_by(Aneuploidy, Dispersal) %>%
#   summarise(pct_pgdis_match = mean(f_pgdis_match)*100,
#             sd_pgdis_match = sd(f_pgdis_match)*100)
# # 
# # 
# ggplot(aneuploidy.values, aes(x = Aneuploidy*100, y = pct_pgdis_match)) +
#   geom_col() +
#   geom_errorbar(ymin = aneuploidy.values$pct_pgdis_match - aneuploidy.values$sd_pgdis_match,
#                 ymax = aneuploidy.values$pct_pgdis_match + aneuploidy.values$sd_pgdis_match,
#                 width = 0.6) +
#   geom_vline(xintercept = 19.5, col="black")+
#   geom_vline(xintercept = 39.5, col="black")+
#   geom_vline(xintercept = 80.5, col="black")+
#   labs(x = "Aneuploidy",
#        y = "Percent of biopsies in correct\nPGDIS class") +
#   coord_cartesian(ylim = c(0, 100))+
#   scale_x_continuous(breaks = seq(0, 100, 10)) +
#   facet_wrap(~Dispersal)+
#   theme_classic()


# dispersal.values = agg.values %>% 
#   filter(Aneuploidy %in% c(0.1, 0.25, 0.5, 0.75, 0.9),
#          Biopsy_size==5) %>%
#   group_by(Aneuploidy, Dispersal) %>%
#   summarise(pct_pgdis_match = mean(f_pgdis_match)*100,
#             pct_merge_match = mean(f_merge_match)*100,
#             sd_pgdis_match = sd(f_pgdis_match)*100,
#             sd_merge_match = sd(f_merge_match)*100)
# 
# ggplot(dispersal.values, aes(x = Dispersal*100, y = pct_merge_match)) +
#   geom_col() + 
#   geom_errorbar(ymin = dispersal.values$pct_pgdis_match - dispersal.values$sd_pgdis_match, 
#                 ymax = dispersal.values$pct_pgdis_match + dispersal.values$sd_pgdis_match,
#                 width = 0.6) +
#   labs(x = "Dispersal of aneuploid cells",
#        y = "Percent of biopsies in correct\nPGDIS class") + 
#   coord_cartesian(ylim = c(0, 100))+
#   scale_x_continuous(breaks = seq(0, 100, 10)) +
#   facet_wrap(~Aneuploidy)+
#   theme_classic()

aneuploidy.dispersal.values = agg_values_5 %>% 
  group_by(Aneuploidy, Dispersal, Biopsy_size) %>%
  summarise(pct_pgdis_match = mean(f_pgdis_match)*100,
            sd_pgdis_match = sd(f_pgdis_match)*100)

aneuploidy.dispersal.heatmap = ggplot(aneuploidy.dispersal.values, aes(x = Aneuploidy*100, y = Dispersal, fill=pct_pgdis_match)) +
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
       fill = "Percent of biopsies\nin correct\nPGDIS class") + 
  coord_cartesian(ylim = c(0, 1.05))+
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_fill_viridis_c()+
  theme_classic()
saveRDS(aneuploidy.dispersal.heatmap, "data/aneuploidy.dispersal.heatmap.Rds")

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



