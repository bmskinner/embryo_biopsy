# Common functions
library(tessera)
library(tidyverse)
library(svglite)
library(data.table)
library(parallel)

# Return the PGDIS class for a given fractional aneuploidy
to.pgdis.class = function(f.aneuploidy){
  case_when(f.aneuploidy < 0.20 ~ "Euploid",
            f.aneuploidy < 0.40 ~ "Low level",
            f.aneuploidy <= 0.80 ~ "High level",
            f.aneuploidy <= 1 ~ "Aneuploid")
}

# Return the merged class for a given fractional aneuploidy
to.merged.class = function(f.aneuploidy){
  case_when(f.aneuploidy < 0.20 ~ "Euploid",
            f.aneuploidy <= 0.80 ~ "Mosaic",
            f.aneuploidy <= 1 ~ "Aneuploid")
}

# Return the mergec class for a given fractional aneuploidy
to.merged.class.2 = function(f.aneuploidy){
  case_when(f.aneuploidy < 0.10 ~ "Euploid",
            f.aneuploidy <= 0.90 ~ "Mosaic",
            f.aneuploidy <= 1 ~ "Aneuploid")
}

# Save a plot as single column width (85mm)
save.single.width = function(plot, filename, height){
  ggsave(plot, filename=paste0(filename, ".png"), dpi=300, units = "mm", width = 85, height = height)
  ggsave(plot, filename=paste0(filename, ".svg"), dpi=300, units = "mm", width = 85, height = height)
}

# Save a plot as double column width (170mm)
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

# Plot PGDIS/merged classification accuracy column charts
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
