# Test ideas for making a random distribution of aneuploid cells
# library(tidyverse)
source("fibonacci.R")
library(plotly)

d = create.blastocyst(n.cells = 200, prop.aneuploid = 0.2, dispersion = 0.5)

plot_ly(x=d$x, y=d$y, z=d$z, 
        type="scatter3d",
        mode="markers",
        color=d$isSeed,
        colors = c("#00FF00", "#FF0000")) %>% layout(showlegend = FALSE)
