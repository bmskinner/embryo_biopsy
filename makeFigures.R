# Make other figures for the paper not created in other analysis files
library(tessera)
library(tidyverse)
library(patchwork)
library(svglite)

save.single.width = function(plot, filename, height){
  ggsave(plot, filename=filename, dpi=300, units = "mm", width = 85, height = height)
}
save.double.width = function(plot, filename, height){
  ggsave(plot, filename=filename, dpi=300, units = "mm", width = 170, height = height)
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
save.single.width(p, filename="figure/Figure_1CD.svg", height = 40)
