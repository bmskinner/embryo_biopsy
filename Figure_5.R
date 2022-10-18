# Figure 5: likely biopsy origins
library(tessera)
library(tidyverse)
library(patchwork)
library(svglite)
library(data.table)
library(parallel)

source("parameters.R")
source("functions.R")

biopsy.values <- do.call(rbind, mclapply(list.files(
  path = AGGREGATE.DATA.PATH,
  pattern = "biopsy", full.names = T
),
fread,
header = T, mc.cores = N.CORES
))

# Calculate the proportion of biopsies originating from each embryo combination
biopsy.aggregate <- biopsy.values %>%
  dplyr::filter(Biopsy_size == 5) %>%
  dplyr::select(
    -embryo_merge_class, -biopsy_merge_class, -n_biopsies_merge,
    -embryo_pgdis_class, -biopsy_pgdis_class
  ) %>%
  dplyr::group_by(Biopsy_size, n_aneuploid) %>%
  dplyr::mutate(
    total_biopsies = sum(n_biopsies_pgdis),
    pct_biopsies = n_biopsies_pgdis / total_biopsies * 100
  )

zero.data <- expand.grid(
  Aneuploidy = ANEUPLOIDY.RANGE,
  Dispersal = DISPERSAL.RANGE,
  pct_biopsies = 0
)

p <- ggplot(biopsy.aggregate, aes(x = Aneuploidy*100, y = Dispersal, fill = pct_biopsies)) +
  geom_raster(data = zero.data) +
  geom_raster() +
  labs(fill = "Percent of\nbiopsies", x = "Embryo aneuploidy (%)", y = "Embryo dispersal") +
  scale_fill_viridis_c() +
  facet_wrap(~n_aneuploid) +
  theme_classic()

save.double.width(p, filename = paste0(FIGURE.OUTPUT.DIR, "/Figure_5_origins"), height = 120)
