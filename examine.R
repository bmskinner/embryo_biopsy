# examine output data

data = read_tsv("conditions.tsv")
filt = data %>% dplyr::rowwise() %>%
  dplyr::mutate(output = str_split(output, ", ")) %>%
  dplyr::mutate(n.euploid   = length(output[output=="0"])) %>%
  dplyr::mutate(n.aneuploid = length(output)-n.euploid,
                PctEuploid  = n.euploid/length(output)*100
                ) %>%
  dplyr::group_by(disps, aneuploids, cells, samples) %>%
  dplyr::mutate(MeanEuploid = mean(PctEuploid),
                SdEuploid   = sd(PctEuploid))


ggplot(filt, aes(x=samples, y=aneuploids, fill=MeanEuploid))+
  geom_tile()+
  scale_fill_viridis_c()+
  facet_grid(disps~cells)
