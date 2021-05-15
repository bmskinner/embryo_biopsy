# Test ideas for making a random disribution of aneuploid cells
library(tidyverse)

calc.sample = function(seed, prop.aneuploid, n.cells){
  n.aneuploid = floor(n.cells*prop.aneuploid)
  
  # The function as used in the excel spreadsheet
  make.vals = function(total, seed, a, c, m){
    vals = list()
    n=seed
    for(i in 1:total){
      colC = (n*a + c)/m
      colD = floor(colC)
      n = (colC-colD) * m
      vals =  c(vals, colC-colD)
    }
    return(unlist(vals))
  }
  
  # Make coordinates and assign type
  data = data.frame("x" = make.vals(n.cells, seed, 785, 3658, 9999), 
                    "y" = make.vals(n.cells, seed, 106, 1283, 9999),
                    "t" = c(rep("a", n.aneuploid), rep("e", n.cells-n.aneuploid))
                    )

  # Take all cells in the rectangular region
  
  
  filt = data %>% filter(x>=0.4 & x<= 0.6 & y>= 0.4 & y<= 0.6)%>% 
    head(n=5)
  
  euploid = filt %>% filter(t=="e") %>% nrow
    
  # ggplot(data, aes(x=x, y=y, col=t))+
  #   geom_point(size=2)+
  #   labs(title=paste0("In region: ", euploid, " Seed: ", seed, " Prop: ", prop.aneuploid, " Cells: ", n.cells))+
  #   theme_classic()+
  #   ggsave(paste0(seed, "_", prop.aneuploid, "_", n.cells, ".png"))
  
  return(euploid)
}

seeds = 1:100
props = seq(0, 1, 0.1)
# cells = seq(200, 1000, by=100)
combinations = expand.grid(seeds, props)

result = data.frame("seed"=NA, "prop"=NA, "count"=NA)
for(r in 1:nrow(combinations)){
  z = calc.sample(combinations$Var1[r], combinations$Var2[r], n.cells=200)
  row.data = c(combinations$Var1[r], combinations$Var2[r], z)
  result = rbind(result, row.data)
}
result  = na.omit(result)

ggplot(result, aes(x=count))+
  geom_bar()+
  facet_wrap(~prop)+
  labs(x="Number of euploid cells sampled", 
       y="Simulations in which this occurred",
       title="Fraction of aneuploid cells in simulations")


ggsave("Proportions.png")
