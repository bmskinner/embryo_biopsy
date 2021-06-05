# examine output data
library(tidyverse)

data = read_tsv("conditions.tsv")
filt = data %>% dplyr::rowwise() %>%
  dplyr::mutate(output = str_split(output, ", ")) %>%
  dplyr::mutate(n.euploid   = length(output[output=="0"])) %>%
  dplyr::mutate(n.aneuploid = length(output)-n.euploid,
                PctEuploid  = n.euploid/length(output)*100
  ) %>%
  dplyr::group_by(disps, aneuploids, cells, samples) %>%
  dplyr::summarise(MeanEuploid = mean(PctEuploid),
                   SdEuploid   = sd(PctEuploid)) 

# ggplot(filt, aes(x=aneuploids, y=MeanEuploid, fill=as.factor(disps), col=as.factor(disps)))+
#   geom_ribbon(aes(ymin=0, ymax=MeanEuploid))+
#   labs(x="Proportion of aneuploid cells", y="All euploid biopsies")+
#   scale_color_viridis_d()+
#   scale_fill_viridis_d()+
#   theme_classic()+
#   facet_grid(samples~cells)


ggplot(filt %>% filter(aneuploids>0 & aneuploids<1), aes(x=cells, y=samples, fill=MeanEuploid>=50))+
  geom_tile()+
  scale_fill_viridis_d()+
  labs(x="Cells in blastocyst", y="Cells biopsied", fill=">=50% chance of finding all euploid")+
  facet_grid(disps~aneuploids)
# 
# 
# ggplot(filt, aes(x=aneuploids, y=MeanEuploid, col=samples, group=aneuploids))+
#   geom_point()+
#   geom_line()+
#   facet_grid(disps~cells)
