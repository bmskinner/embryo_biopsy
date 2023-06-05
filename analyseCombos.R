# Analyse pre-computed raw biopsy data to create aggregate values for charting
library(data.table)
library(tidyverse)
library(parallel)
library(fs)

source("parameters.R")
source("functions.R")

# Given raw input data, calculate the accuracy of the biopsies
# with respect to embryo classification boundaries
calc.pgdis.accuracy <- function(data) {
  data %>%
    rowwise() %>%
    mutate(
      n_aneuploid = list(c_across(starts_with("V"))),
      f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size))
    ) %>%
    select(-starts_with("V")) %>%
    mutate(
      pgdis_class = list(to.pgdis.class(f_aneuploid)),
      equal_class = list(to.equal.class(f_aneuploid)),
      merge_class = list(to.merged.class(f_aneuploid)),
      actual_pgdis_class = to.pgdis.class(Aneuploidy),
      actual_equal_class = to.equal.class(Aneuploidy),
      actual_merge_class = to.merged.class(Aneuploidy),
      n_pgdis_match = sum(pgdis_class == actual_pgdis_class),
      f_pgdis_match = n_pgdis_match / 200,
      n_equal_match = sum(equal_class == actual_equal_class),
      f_equal_match = n_equal_match / 200,
      n_merge_match = sum(merge_class == actual_merge_class),
      f_merge_match = n_merge_match / 200
    ) %>%
    select(-n_aneuploid, -f_aneuploid, -pgdis_class, -merge_class, -equal_class) %>%
    ungroup()
}

# Calculate the number of biopsies in each class for each aneuploidy, dispersal
# and biopsy size combination
calc.biopsy.accuracy <- function(data) {
  data %>%
    rowwise() %>%
    dplyr::mutate(n_aneuploid = list(c_across(starts_with("V")))) %>%
    dplyr::select(-starts_with("V")) %>%
    tidyr::unnest(n_aneuploid) %>%
    dplyr::mutate(
      embryo_pgdis_class = to.pgdis.class(Aneuploidy),
      embryo_equal_class = to.equal.class(Aneuploidy),
      embryo_merge_class = to.merged.class(Aneuploidy),
      biopsy_pgdis_class = to.pgdis.class(n_aneuploid / Biopsy_size),
      biopsy_equal_class = to.equal.class(n_aneuploid / Biopsy_size),
      biopsy_merge_class = to.merged.class(n_aneuploid / Biopsy_size)
    ) %>%
    dplyr::ungroup() %>% # cancels rowwise
    dplyr::group_by(Aneuploidy, Dispersal, n_aneuploid, embryo_pgdis_class, biopsy_pgdis_class, Biopsy_size) %>%
    dplyr::mutate(n_biopsies_pgdis = dplyr::n()) %>%
    dplyr::ungroup() %>% # cancels rowwise
    dplyr::group_by(Aneuploidy, Dispersal, n_aneuploid, embryo_equal_class, biopsy_equal_class, Biopsy_size) %>%
    dplyr::mutate(n_biopsies_equal = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Aneuploidy, Dispersal, n_aneuploid, embryo_merge_class, biopsy_merge_class, Biopsy_size) %>%
    dplyr::mutate(n_biopsies_merge = dplyr::n()) %>%
    dplyr::select(-Seed) %>%
    dplyr::distinct()
}

# Calculate aggregate data from raw values for making heatmaps
make.aggregate.values <- function(e, a, d) {
  aneu.part.file <- paste0(AGGREGATE.DATA.PATH, "/merged_e", e, "_a", a, "_d", d, ".csv")
  if (!file.exists(aneu.part.file)) {
    in.file <- paste0(RAW.DATA.PATH, "/raw_values_e", e, "_a", a, "_d", d, ".csv")
    if (file.exists(in.file)) {
      in.data <- fread(file = in.file, header = T)
      filt.tf <- calc.pgdis.accuracy(in.data)
      write.csv(filt.tf, file = aneu.part.file, quote = F, row.names = F)
      rm(filt.tf)
      gc()
    }
  }
}

# Make output files aggregating data for biopsy accuracy calculations
make.biopsy.values <- function(e, a) {
  out.file <- paste0(AGGREGATE.DATA.PATH, "/biopsy_e", e, "_a", a, ".csv")
  if (!file.exists(out.file)) {
    in.files <- list.files(
      path = RAW.DATA.PATH,
      pattern = paste0("raw_values_e", e, "_a", a, "_"), full.names = T
    )

    in.data <- do.call(rbind, mclapply(in.files,
      fread,
      header = T,
      mc.cores = N.CORES
    ))
    filt.tf <- calc.biopsy.accuracy(in.data)
    write.csv(filt.tf, file = out.file, quote = F, row.names = F)
    rm(filt.tf)
    invisible(gc(verbose = F))
  }
}


# Make output directories
if (!fs::dir_exists(AGGREGATE.DATA.PATH)) {
  fs::dir_create(AGGREGATE.DATA.PATH, recursive = T)
}
if (!fs::dir_exists(RAW.DATA.PATH)) {
  fs::dir_create(RAW.DATA.PATH, recursive = T)
}

# Create the aggregate files
combinations <- expand.grid(a = ANEUPLOIDY.RANGE, d = DISPERSAL.RANGE, e = EMBRYO.SIZES)
mcmapply(make.aggregate.values, e = combinations$e, a = combinations$a, d = combinations$d, mc.cores = N.CORES)

combinations <- expand.grid(a = ANEUPLOIDY.RANGE, e = EMBRYO.SIZES)
mcmapply(make.biopsy.values, e = combinations$e, a = combinations$a, mc.cores = N.CORES)
