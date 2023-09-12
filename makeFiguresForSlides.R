# Generate charts for slides
library(tidyverse)

aggregate.data <- do.call(rbind, lapply(
  list.files(path = "data/aggregates", pattern = "*.csv", full.names = T),
  function(f) read.csv(f, header = T)
))


biopsy.5.data <- aggregate.data[aggregate.data$Biopsy_size == 5, ]


dispersal.aneuploidy.heatmap <- ggplot(biopsy.5.data, aes(
  x = Aneuploidy * 100,
  y = Dispersal,
  fill = f_pgdis_match * 100
)) +
  geom_raster() +
  scale_fill_viridis_c(breaks = seq(0, 100, 20)) +
  labs(
    x = "Percent aneuploidy in embryo",
    y = "Dispersal of aneuploid cells",
    fill = "Percentage\nof accurate\nPGDIS\nclassifications"
  ) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  annotate(geom = "text", label = "Euploid", x = 10, y = 1.06) +
  annotate(geom = "text", label = "Low", x = 30, y = 1.06) +
  annotate(geom = "text", label = "High", x = 60, y = 1.06) +
  annotate(geom = "text", label = "Aneuploid", x = 90, y = 1.06) +
  geom_vline(xintercept = 19.5, col = "white") +
  geom_vline(xintercept = 39.5, col = "white") +
  geom_vline(xintercept = 80.5, col = "white") +
  theme_classic()
saveRDS(dispersal.aneuploidy.heatmap, "data/dispersal.aneuploidy.heatmap.Rds")


# dispersal.aneuploidy.biopsy.heatmap = ggplot(aggregate.data[aggregate.data$Biopsy_size %in% c(3, 5, 10, 20),], aes(x = Aneuploidy*100,
#                                                          y = Dispersal,
#                                                          fill = f_pgdis_match*100)) +
#   geom_raster() +
#   scale_fill_viridis_c(breaks = seq(0, 100, 20)) +
#   labs(x = "Percent aneuploidy in embryo",
#        y = "Dispersal of aneuploid cells",
#        fill = "Percentage\nof accurate\nPGDIS\nclassifications") +
#   scale_x_continuous(breaks = seq(0, 100, 10)) +
#   scale_y_continuous(breaks = seq(0, 1, 0.1)) +
#   annotate(geom = "text", label = "Euploid", x=10, y=1.06)+
#   annotate(geom = "text", label = "Low", x=30, y=1.06)+
#   annotate(geom = "text", label = "High", x=60, y=1.06)+
#   annotate(geom = "text", label = "Aneuploid", x=90, y=1.06)+
#   geom_vline(xintercept = 19.5, col="white")+
#   geom_vline(xintercept = 39.5, col="white")+
#   geom_vline(xintercept = 80.5, col="white")+
#   facet_wrap(~Biopsy_size, ncol = 2)+
#   theme_classic()
# saveRDS(dispersal.aneuploidy.biopsy.heatmap, "data/dispersal.aneuploidy.biopsy.heatmap.Rds")
