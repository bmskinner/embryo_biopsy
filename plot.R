source("grid.R")
library(tidyverse)
####################################
# Running the analysis
####################################

runAnalysis = function(){
  # Make a grid of parameters
  props = seq(0.0, 1.0, 0.01)
  disps = seq(0.0, 1.0, 0.05) # dispersion - if 0, all one clump, if 1 all cells individual
  combs = expand.grid(props, disps)
  
  # Prepare result data frame
  result = data.frame("prop"=NA, "disp"=NA, "e"=NA)
  
  iterations = 500 # number of simulations per parameter combination
  
  for(i in 1:nrow(combs)){
    n.aneu = floor(combs$Var1[i] * 10 * 10)
    cat("Combination", i, "of", nrow(combs), "Disp:", combs$Var2[i], "Prop:", combs$Var1[i],"for",n.aneu, "aneuploid cells\n")
    ptm <- proc.time()
    counts.and.plots = lapply(1:iterations, make.count,
                              prop.aneuploids=combs$Var1[i],
                              dispersion=combs$Var2[i], make.chart=F, dim.x=10, dim.y=10)
    time.taken = proc.time() - ptm
    cat("Combination", i, "took ", time.taken[['elapsed']], "seconds\n")
    # Extract counts and plots separately
    counts = sapply(1:iterations, function(i) counts.and.plots[[i]][['n.euploid']])
    plots  = lapply(1:iterations, function(i) counts.and.plots[[i]][['plot']])
    
    temp = data.frame("prop"=rep(combs$Var1[i], iterations),
                      "disp"=rep(combs$Var2[i], iterations),
                      "e" = counts)
    result = rbind(result, temp)
  }
  
  # Plot the results
  result = na.omit(result)
  
  bar.plot = ggplot(result, aes(x=e))+
    geom_bar()+
    facet_grid(prop~disp)+
    labs(title=paste("Cols = dispersal (0=clumped, 1=dispersed), rows = proportion aneuploid cells"),
         x="Number of euploid cells in 5 cell biopsy",
         y="Simulations in which this occurred")
  
  
  vals = result %>% group_by(prop, disp) %>%
    summarise(meanEuploid = mean(e))
  
  return(list('bars' = bar.plot, 'vals'=vals))
  
}

result = runAnalysis()

ggplot(result[['vals']], aes(x=disp, y=prop, fill=meanEuploid))+
  geom_tile()+
  labs(x="Dispersal", y="Proportion aneuploid", fill="Mean euploid")+
  scale_fill_viridis_c()+
  theme_minimal()
ggsave("Prop_by_dispersal.png")