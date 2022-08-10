# Check the values from Ren et al
library(tidyverse)
ren.data = read.csv("Ren/Supp_tables_merged.csv", header = T)

ren.data$isEuploid = ren.data$Karyotype=="46,XY" | ren.data$Karyotype=="46,XX"

ren.filt = ren.data %>% dplyr::group_by(Embryoid, Type) %>%
  dplyr::mutate(TotalCells = dplyr::n()) %>%
  dplyr::group_by(Embryoid, Type, Karyotype) %>%
  dplyr::mutate(Count = dplyr::n(), PctTotal = Count/TotalCells) %>%
  dplyr::select(Embryoid, Type, Karyotype, Count, TotalCells, PctTotal) %>%
  dplyr::distinct() %>%
  dplyr::filter(Type == "TE")
  

cat("Total embryos:", length(unique(ren.filt$Embryoid)))


ren.wide = ren.filt %>% 
  dplyr::select(-PctTotal) %>% 
  tidyr::pivot_wider(names_from = Karyotype, values_from = Count, values_fill = 0)

#  From TE biopsy alone:
n.euploid = nrow(ren.wide[ren.wide$`FALSE`==0,])
n.aneuploid = nrow(ren.wide[ren.wide$`TRUE`==0,])
n.mosaic = nrow(ren.wide) - n.euploid -n.aneuploid
