# Make other figures for the paper not created in other analysis files
library(tessera)
library(tidyverse)
library(patchwork)
library(svglite)

# source("analyseCombos.R")

save.single.width = function(plot, filename, height){
  ggsave(plot, filename=paste0(filename, ".png"), dpi=300, units = "mm", width = 85, height = height)
  ggsave(plot, filename=paste0(filename, ".svg"), dpi=300, units = "mm", width = 85, height = height)
}
save.double.width = function(plot, filename, height){
  ggsave(plot, filename=paste0(filename, ".png"), dpi=300, units = "mm", width = 170, height = height)
  ggsave(plot, filename=paste0(filename, ".svg"), dpi=300, units = "mm", width = 170, height = height)
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




