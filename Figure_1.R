# Figure 1: embryo models and their biopsy distribution

library(tessera)
library(tidyverse)
library(patchwork)
library(svglite)
library(data.table)

source("parameters.R")
source("functions.R")


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
save.single.width(p, filename=paste0(FIGURE.OUTPUT.DIR, "/Figure_1CD"), height = 40)