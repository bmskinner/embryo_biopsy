# Check the values from Ren et al
library(tidyverse)
ren.data = read.csv("Ren/Supp_tables_merged.csv", header = T)

ren.data$isEuploid = ren.data$Karyotype=="46,XY" | ren.data$Karyotype=="46,XX"

# Process to match Ren et al's fig S3.
# Note that because there are ties in the number of cells of different karyotypes,
# we arbirarily select the karyotype to be 'primary' when calculating mosaicism rates.
# This is why we use rownumber and not max(Count)
ren.filt = ren.data %>% dplyr::group_by(Embryoid, Type) %>%
  dplyr::mutate(TotalCells = dplyr::n()) %>%
  dplyr::group_by(Embryoid, Type, Karyotype, isEuploid) %>%
  dplyr::mutate(Count = dplyr::n(), PctTotal = Count/TotalCells) %>%
  dplyr::select(Embryoid, Type, Karyotype, Count, TotalCells, PctTotal, isEuploid) %>%
  dplyr::distinct() %>%
  dplyr::group_by(Embryoid, Type) %>%
  dplyr::arrange(Embryoid, Type, Count, by_group=T ) %>% # order so the highest count is last
  dplyr::mutate(isMosaic = length(unique(Karyotype))>1,
                isPrimaryKaryotype = row_number() == max(row_number())) %>% #set the last row of the group to be the primary karyotype
  dplyr::filter(Type == "TE")%>% 
  dplyr::ungroup() %>%
  dplyr::group_by(Embryoid, Type) %>%
  dplyr::mutate(anyEuploid = any(isEuploid),
                allEuploid = all(isEuploid))
  

#  From TE biopsy alone:
n.total = length(unique(ren.filt$Embryoid))
n.euploid = length(unique(ren.filt[ren.filt$isMosaic==F & ren.filt$isEuploid==T,]$Embryoid))
n.aneuploid = length(unique(ren.filt[ren.filt$isMosaic==F & ren.filt$isEuploid==F,]$Embryoid))
n.mosaic = length(unique(ren.filt[ren.filt$isMosaic==T,]$Embryoid))

cat("Total embryos    :", n.total, "\nMosaic embryos   :", n.mosaic, "\nEuploid embryos  :", n.euploid,"\nAneuploid embryos:", n.aneuploid, "\n")

ren.wide = ren.filt %>% 
  dplyr::filter(isMosaic==T & isPrimaryKaryotype==F) %>%
  dplyr::group_by(Embryoid) %>%
  dplyr::summarise(MosaicismRate = sum(PctTotal), n.cells = unique(TotalCells))

# Digitised copy of Fig S3 TE section
ren.transcribed = read.csv("Ren/ren.te.data.csv", header = T) %>% dplyr::mutate(n_cells_sequenced = round(n_cells_sequenced))

# Compare our calculations to Rens
ggplot(ren.wide, aes(x = n.cells, y=MosaicismRate))+
  geom_point(size = 4)+
  geom_point(data = ren.transcribed, aes(x=n_cells_sequenced, y=mosaicism_rate_pct/100), col="red")+
  scale_x_continuous(breaks = seq(2, 15, 2))+
  scale_y_continuous(breaks = seq(0, 1, 0.1))
# Should have a good overlap at this point
                     

# Note that this is measuring the level of mosaicism relative to the primary karyotype, not to euploidy.
# We need to recode these for our own analysis.

our.euploid = ren.filt[ren.filt$isEuploid==T & ren.filt$isMosaic==F,]

our.mosaic = ren.filt %>%
  dplyr::filter(anyEuploid & !allEuploid & !isEuploid) %>%
  dplyr::summarise(PctAenuploid = sum(Count)/TotalCells) %>%
  dplyr::distinct()

our.aneuploid = ren.filt %>%
  dplyr::filter(!anyEuploid)  %>%
  dplyr::summarise(PctAenuploid = sum(Count)/TotalCells) %>%
  dplyr::distinct()

cat("Total embryos    :", n.total, "\nMosaic embryos   :", nrow(our.mosaic), "\nEuploid embryos  :", nrow(our.euploid),"\nAneuploid embryos:", nrow(our.aneuploid), "\n")
# We can now use the PctAneuploid values from the mosaic embryos in our models
