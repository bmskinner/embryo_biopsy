# Create two cell biopsies

library(tessera)
library(parallel)
library(tidyverse)
library(data.table)

ANEUPLOIDY.RANGE = seq(0, 1, 0.01)
DISPERSAL.RANGE = seq(0, 1, 0.5)
N.REPLICATES = 100 
BIOPSY.SIZES = c(3:10, 15, 20, 25, 30)
EMBRYO.SIZES = seq(100, 250, 50)
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

# For a given embryo, take all possible two biopsies
take.two.biopsies = function(embryo, biopsy.size, chromosome){
  
  # cat("Taking two biopsies\n")
  getNeighboursToExclude = function(embryo, cell.index){
    # get the direct neighbours - these will be in the biopsy
    neighbours = tessera::getNeighbouringCellIndexes(embryo, cell.index)
    
    # Get neighbours of neighbours so we don't have any surrounding biopsies 
    # including the missing cells
    neighbours2 = sapply(neighbours, tessera::getNeighbouringCellIndexes, embryo=embryo)
    return(unique(c(cell.index, neighbours, neighbours2)))
  }
  
  all.cells = 1:length(embryo)
  
  # Get all biopsies in a vector indexed by cell
  all.biopsies = sapply(all.cells, function(i){
    tessera::takeBiopsy(embryo, biopsy.size = biopsy.size, chromosome = chromosome, index.cell = i)
  })
  
  # Now find the biopsy combinations to look up
  unlist(sapply(1:length(embryo), function(i){
    # We now need to exclude neighbouring cells for second biopsy
    # Cells that are not neighbours or neighbours of neighbours
    keep.list = all.cells[!all.cells %in% getNeighboursToExclude(embryo, i)]
    
    sapply(keep.list, function(j){ (all.biopsies[i] + all.biopsies[j])/2 })
  }))
}

# For a combination of aneuploidy and dispersal, calculate all two biopsy values
make.two.biopsy.data = function(aneuploidy, dispersal){
  
  part.file = paste0("data/two_biopsy/two_biopsy_values_a", aneuploidy, "_d", dispersal, ".csv")
  
  if(!file.exists(part.file)){
    
    get.biopsies = function(s, embryo.size, biopsy.size){
      e = tessera::Embryo(n.cells = embryo.size,
                          prop.aneuploid = aneuploidy,
                          dispersal = dispersal,
                          concordance = 1, rng.seed = s)
      b = take.two.biopsies(e, biopsy.size = biopsy.size, chromosome = 1)
      
      list(Aneuploidy = aneuploidy, Dispersal = dispersal, 
           Embryo_size = embryo.size, Biopsy_size = biopsy.size, n_aneuploid = b)
    }
    
    ad.combo = expand.grid(embryo.size  = EMBRYO.SIZES,
                           biopsy.size  = BIOPSY.SIZES,
                           seed         = 1:N.REPLICATES)
    
    r = do.call(rbind, mcmapply(get.biopsies, 
                                s           = ad.combo$seed, 
                                embryo.size = ad.combo$embryo.size,
                                biopsy.size = ad.combo$biopsy.size,
                                mc.cores = N.CORES, SIMPLIFY = F)) %>% 
      as.data.frame %>%
      tidyr::unnest_longer(., col=n_aneuploid, simplify = T)
      
    r$Aneuploidy  = unlist(r$Aneuploidy)
    r$Dispersal   = unlist(r$Dispersal)
    r$Embryo_size = unlist(r$Embryo_size)
    r$Biopsy_size = unlist(r$Biopsy_size)
    
    r = r %>% 
      dplyr::group_by(Aneuploidy, Dispersal, Embryo_size, Biopsy_size, n_aneuploid) %>%
      dplyr::summarise(Total_biopsies = dplyr::n())
    
    write.table(r, part.file, append = F, col.names = T, row.names = F, sep = "\t", quote = F)
    rm(r)
    gc()
  }
}


combinations = expand.grid(aneuploidies = ANEUPLOIDY.RANGE,
                           dispersals   = DISPERSAL.RANGE)

# Make the combinations
mapply(make.two.biopsy.data, 
       aneuploidy = combinations$aneuploidies,
       dispersal  = combinations$dispersals)


# Read and analyse the saved data
in.files = list.files(path = "data/two_biopsy",
                      pattern = paste0("two_biopsy"), full.names = T)

in.data = do.call(rbind, mclapply(in.files,
                                  fread,
                                  header = T,
                                  mc.cores = 22))

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
  ggplot(filt.data, aes(x = f_aneuploid, y = Aneuploidy, fill=PctTotal))+
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

  ggsave(filename = paste0("figure/predictive_heatmap_two_biopsy_pgdis_",b,".png"), 
         dpi=300, units = "mm", width = 170, height = 150)
  
  # Plot with merge classes
  ggplot(filt.data, aes(x = f_aneuploid, y = Aneuploidy, fill=PctTotal))+
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
  ggsave(filename = paste0("figure/predictive_heatmap_two_biopsy_merge_",b,".png"), 
         dpi=300, units = "mm", width = 170, height = 150)

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
  
  ggplot(col.data.pgdis, aes(x = n_aneuploid/b*100, y= PctTotal,  fill=IsPGDISCorrect))+
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
  ggsave(filename = paste0("figure/predictive_columns_two_biopsy_pgdis_",b,".png"), 
         dpi=300, units = "mm", width = 170, height = 150)
  
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
  
  ggplot(col.data.merge, aes(x = n_aneuploid/b*100, y= PctTotal,  fill=IsPGDISCorrect))+
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
  ggsave(filename = paste0("figure/predictive_columns_two_biopsy_merge_",b,".png"),
         dpi=300, units = "mm", width = 170, height = 150)
  
}
