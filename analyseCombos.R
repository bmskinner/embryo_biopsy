# Test analysing pre-computed values
library(data.table)
library(tidyverse)
library(parallel)

ANEUPLOIDY.RANGE = seq(0, 1, 0.01)
DISPERSAL.RANGE = seq(0, 1, 0.01)
N.REPLICATES = 100 
BIOPSY.SIZES = c(3:10, 15, 20, 25, 30)
N.CORES = ifelse(Sys.info()["sysname"]=="Windows", 1, 20) 

data.file = "data/all.combos.csv"

# Given raw input data, calculate the difference between the biopsy and 
# embryo aneuploidies and summarise into counts
calc.step.accuracy = function(data, var.name){
  data %>%
    rowwise %>%
    mutate(n_aneuploid = list(c_across(starts_with("V"))),
           f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size)),
           diff_embryo = list(abs(f_aneuploid-Aneuploidy))) %>%
    select(-starts_with("V"), -n_aneuploid, -f_aneuploid) %>%
    ungroup %>%
    unnest(diff_embryo) %>%
    dplyr::group_by(Aneuploidy, Dispersal, Biopsy_size, diff_embryo) %>%
    dplyr::summarise(Count = dplyr::n())
}

# Given raw input data, calculate the accuracy of the biopsies
# with respect to PGDIS classification boundaries
calc.pgdis.accuracy = function(data){
  
  to.pgdis.class = function(f.aneuploidy){
    case_when(f.aneuploidy < 0.20 ~ "Euploid",
              f.aneuploidy < 0.40 ~ "Low level",
              f.aneuploidy <= 0.80 ~ "High level",
              f.aneuploidy <= 1 ~ "Aneuploid")
  }
  
  to.merged.class = function(f.aneuploidy){
    case_when(f.aneuploidy < 0.20 ~ "Euploid",
              f.aneuploidy <= 0.80 ~ "Mosaic",
              f.aneuploidy <= 1 ~ "Aneuploid")
  }
  
  data %>%
    rowwise %>%
    mutate(n_aneuploid = list(c_across(starts_with("V"))),
           f_aneuploid = list(as.double(c_across(starts_with("V")) / Biopsy_size))) %>%
    select(-starts_with("V")) %>%
    mutate(pgdis_class = list(to.pgdis.class(f_aneuploid)),
           merge_class = list(to.merged.class(f_aneuploid)),
           actual_pgdis_class = to.pgdis.class(Aneuploidy),
           actual_merge_class = to.merged.class(Aneuploidy),
           n_pgdis_match = sum(pgdis_class == actual_pgdis_class),
           f_pgdis_match = n_pgdis_match / 200,
           n_merge_match = sum(merge_class==actual_merge_class),
           f_merge_match = n_pgdis_match / 200) %>%
    select(-n_aneuploid, -f_aneuploid, -pgdis_class, -merge_class) %>%
    ungroup
}

# Read the saved raw values
values = data.table::fread(data.file, header = T)


aggregate.values = function(b, a, d){
  aneu.part.file = paste0("data/aggregates/merged_b",b,"_a", a,"_d",d,".csv")
  if(!file.exists(aneu.part.file)){
    cat("Summarising a:", a, "\td: ", d, "\tb:",b,  "\n")
    filt = values[values$Aneuploidy == a & values$Dispersal == d & values$Biopsy_size==b,]
    filt.tf = calc.pgdis.accuracy(filt)
    write.csv(filt.tf, file = aneu.part.file, quote = F, row.names = F)
    rm(values_tf)
    rm(filt)
    gc()
  }
}

combinations = expand.grid(a = ANEUPLOIDY.RANGE, d = DISPERSAL.RANGE, b = BIOPSY.SIZES)

# Function writes output files, no need to store in object
mcmapply(aggregate.values, a=combinations$a, d = combinations$d, b = combinations$b, mc.cores=N.CORES)


# # Summarising is memory intensive, chunk and save
# for(b in BIOPSY.SIZES){
#   
#   for(a in seq(0, 1, 0.25)){
#     e = a+0.25
#     a = ifelse(a==0, -0.1, a) # ensure zero included
#     cat("Summarising ", a, "-", e, "\n")
#     
#     aneu.part.file = paste0("data/aggregates/part_b",b,"_a", a,"-",e,".csv")
#     if(!file.exists(aneu.part.file)){
#       filt = values[values$Aneuploidy > a & values$Aneuploidy <= e & values$Biopsy_size==b,]
#       values.tf = calc.pgdis.accuracy(filt)
#       write.csv(values.tf, file = aneu.part.file, quote = F, row.names = F)
#       rm(values_tf)
#       rm(filt)
#       gc()
#     }
#   }
# }

# Now we can generate the summary figures

agg_values = do.call(rbind, lapply(list.files(path = "data/aggregates", 
                                              pattern = ".csv", full.names = T), 
                                   fread, 
                                   header = T))

agg_values_5 = agg_values[agg_values$Biopsy_size==5,]



aneuploidy.values = agg_values %>% 
  filter(Dispersal %in% c(0, 0.25, 0.5, 0.75, 1),
         Biopsy_size==5) %>%
  group_by(Aneuploidy, Dispersal) %>%
  summarise(pct_pgdis_match = mean(f_pgdis_match)*100,
            sd_pgdis_match = sd(f_pgdis_match)*100)


ggplot(aneuploidy.values, aes(x = Aneuploidy*100, y = pct_pgdis_match)) +
  geom_col() + 
  geom_errorbar(ymin = aneuploidy.values$pct_pgdis_match - aneuploidy.values$sd_pgdis_match, 
                ymax = aneuploidy.values$pct_pgdis_match + aneuploidy.values$sd_pgdis_match,
                width = 0.6) +
  geom_vline(xintercept = 19.5, col="black")+
  geom_vline(xintercept = 39.5, col="black")+
  geom_vline(xintercept = 80.5, col="black")+
  labs(x = "Aneuploidy",
       y = "Percent of biopsies in correct\nPGDIS class") + 
  coord_cartesian(ylim = c(0, 100))+
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  facet_wrap(~Dispersal)+
  theme_classic()


dispersal.values = agg_values %>% 
  filter(Aneuploidy %in% c(0.1, 0.25, 0.5, 0.75, 0.9),
         Biopsy_size==5) %>%
  group_by(Aneuploidy, Dispersal) %>%
  summarise(pct_pgdis_match = mean(f_pgdis_match)*100,
            sd_pgdis_match = sd(f_pgdis_match)*100)

ggplot(dispersal.values, aes(x = Dispersal*100, y = pct_pgdis_match)) +
  geom_col() + 
  geom_errorbar(ymin = dispersal.values$pct_pgdis_match - dispersal.values$sd_pgdis_match, 
                ymax = dispersal.values$pct_pgdis_match + dispersal.values$sd_pgdis_match,
                width = 0.6) +
  labs(x = "Dispersal of aneuploid cells",
       y = "Percent of biopsies in correct\nPGDIS class") + 
  coord_cartesian(ylim = c(0, 100))+
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  facet_wrap(~Aneuploidy)+
  theme_classic()

aneuploidy.dispersal.values = agg_values_5 %>% 
  group_by(Aneuploidy, Dispersal, Biopsy_size) %>%
  summarise(pct_pgdis_match = mean(f_pgdis_match)*100,
            sd_pgdis_match = sd(f_pgdis_match)*100)

aneuploidy.dispersal.heatmap = ggplot(aneuploidy.dispersal.values, aes(x = Aneuploidy*100, y = Dispersal, fill=pct_pgdis_match)) +
  geom_raster()+
  geom_vline(xintercept = 19.5, col="white")+
  geom_vline(xintercept = 39.5, col="white")+
  geom_vline(xintercept = 80.5, col="white")+
  annotate(geom = "text", label = "Euploid", x=10, y=1.06)+
  annotate(geom = "text", label = "Low", x=30, y=1.06)+
  annotate(geom = "text", label = "High", x=60, y=1.06)+
  annotate(geom = "text", label = "Aneuploid", x=90, y=1.06)+
  labs(x = "Aneuploidy",
       y = "Dispersal",
       fill = "Percent of biopsies\nin correct\nPGDIS class") + 
  coord_cartesian(ylim = c(0, 1.05))+
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_fill_viridis_c()+
  theme_classic()
saveRDS(aneuploidy.dispersal.heatmap, "data/aneuploidy.dispersal.heatmap.Rds")

aneuploidy.dispersal.biopsy.values = agg_values %>% 
  group_by(Aneuploidy, Dispersal, Biopsy_size) %>%
  summarise(pct_pgdis_match = mean(f_pgdis_match)*100,
            sd_pgdis_match = sd(f_pgdis_match)*100)

aneuploidy.dispersal.biopsy.heatmap = ggplot(aneuploidy.dispersal.biopsy.values, aes(x = Aneuploidy*100, y = Dispersal, fill=pct_pgdis_match)) +
  geom_raster()+
  geom_vline(xintercept = 19.5, col="white")+
  geom_vline(xintercept = 39.5, col="white")+
  geom_vline(xintercept = 80.5, col="white")+
  labs(x = "Aneuploidy",
       y = "Dispersal",
       fill = "Percent of biopsies\nin correct\nPGDIS class") + 
  coord_cartesian(ylim = c(0, 1))+
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_fill_viridis_c()+
  facet_wrap(~Biopsy_size, ncol = 4)+
  theme_classic()
saveRDS(aneuploidy.dispersal.biopsy.heatmap, "data/aneuploidy.dispersal.biopsy.heatmap.Rds")



