# Analyse pre-computed raw biopsy data to create aggregate values for charting
library(data.table)
library(tidyverse)
library(parallel)
library(svglite)
library(patchwork)

ANEUPLOIDY.RANGE = seq(0, 1, 0.01)
DISPERSAL.RANGE = seq(0, 1, 0.01)
N.REPLICATES = 100 
BIOPSY.SIZES = c(3:10, 15, 20, 25, 30)
EMBRYO.SIZES = c(100, 150, 200, 250)
N.CORES = ifelse(Sys.info()["sysname"]=="Windows", 1, 5) 

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

# Given raw input data, calculate the difference between the biopsy and 
# embryo aneuploidies and summarise into cumulative percentage
calc.step.accuracy = function(data, all.vals){
  
  data %>%
    rowwise %>%
    mutate(f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size))) %>%
    select(-starts_with("V")) %>%
    ungroup %>%
    unnest(f_aneuploid) %>%
    dplyr::mutate(diff_to_embryo = round(abs(f_aneuploid-Aneuploidy), digits=2)) %>%
    dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size, diff_to_embryo) %>%
    dplyr::summarise(Count = dplyr::n()) %>%
    dplyr::bind_rows(., all.vals) %>% # add in all data
    dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size, diff_to_embryo) %>%
    dplyr::summarise(Count = sum(Count)) %>% # recalculate counts
    dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size) %>%
    dplyr::mutate(TotalBiopsies = sum(Count),
                  PctBiopsies   = Count/TotalBiopsies*100) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Biopsy_size, Dispersal, Aneuploidy, diff_to_embryo) %>%
    dplyr::select(-TotalBiopsies, -Count) %>%
    dplyr::distinct() %>% 
    dplyr::mutate(CumPct = cumsum(PctBiopsies))
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
make.aggregate.values = function(e, a, d){
  aneu.part.file = paste0("data/aggregates/merged_e", e, "_a", a, "_d", d, ".csv")
  if(!file.exists(aneu.part.file)){
    in.file = paste0("data/raw/raw_values_e", e, "_a", a, "_d", d, ".csv")
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

combinations = expand.grid(a = ANEUPLOIDY.RANGE, d = DISPERSAL.RANGE, e = EMBRYO.SIZES)

# Functions write output files, no need to store in object
mcmapply(make.aggregate.values, e = combinations$e, a = combinations$a, d = combinations$d, mc.cores = N.CORES)
mcmapply(make.biopsy.values, a = ANEUPLOIDY.RANGE, mc.cores = 3)
