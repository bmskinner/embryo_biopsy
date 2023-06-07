#  Figure 4 - like Fig 3 faceted over all biopsy sizes
library(tidyverse)
library(patchwork)
library(data.table)
library(parallel)

source("parameters.R")
source("functions.R")

# Is the biopsy saying the embryo is better or worse than reality?
find.type <- function(f, a) {
  dplyr::case_when(
    f < a ~ "Better",
    f == a ~ "Equal",
    T ~ "Worse"
  )
}

# Given raw input data, calculate the accuracy of the biopsies
# with respect to embryo classification boundaries
calc.better.worse <- function(data) {
  data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      n_aneuploid = list(c_across(starts_with("V"))),
      f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size))
    ) %>%
    dplyr::select(-starts_with("V")) %>%
    dplyr::mutate(
      type = list(find.type(f_aneuploid, Aneuploidy))
    ) %>%
    tidyr::unnest(type) %>%
    dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size, Embryo_size, type) %>%
    dplyr::mutate(n_type = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size, Embryo_size) %>%
    dplyr::mutate(
      total = dplyr::n(),
      f_type = n_type / total
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-n_aneuploid, -f_aneuploid, -Seed) %>%
    dplyr::distinct()
}

# Calculate aggregate data from raw values for making heatmaps
make.aggregate.values <- function(e, a, d) {
  part.file <- paste0(AGGREGATE.DATA.PATH, "/better.worse_e", e, "_a", a, "_d", d, ".csv")
  in.file <- paste0(RAW.DATA.PATH, "/raw_values_e", e, "_a", a, "_d", d, ".csv")

  if (!file.exists(part.file) & file.exists(in.file)) {
    in.data <- fread(file = in.file, header = T)
    filt.tf <- calc.better.worse(in.data)
    write.csv(filt.tf, file = part.file, quote = F, row.names = F)
    rm(filt.tf)
  }
}

# Create the aggregate files
combinations <- expand.grid(a = ANEUPLOIDY.RANGE, d = DISPERSAL.RANGE, e = EMBRYO.SIZES)

mcmapply(make.aggregate.values, e = combinations$e, a = combinations$a, d = combinations$d, mc.cores = N.CORES)

# Read in the aggregate data and analyse

agg.values <- do.call(rbind, mclapply(list.files(
  path = AGGREGATE.DATA.PATH,
  pattern = "better.worse", full.names = T
),
fread,
header = T, mc.cores = N.CORES
)) %>%
  dplyr::distinct()

p <- ggplot(agg.values[agg.values$Embryo_size == 200, ], aes(x = Aneuploidy, y = Dispersal, fill = f_type)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0, 1)) +
  facet_grid(Biopsy_size ~ type)


save.double.width(p, filename = paste0(FIGURE.OUTPUT.DIR, "/Better_or_worse_heatmap"), height = 170)
