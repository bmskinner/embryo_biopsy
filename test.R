# Test ideas for making a random distribution of aneuploid cells
# library(tidyverse)
source("fibonacci.R")
library(plotly)

d = create.blastocyst(n.cells = 144, prop.aneuploid = 0.25, dispersion = 1)

plot_ly(x=d$x, y=d$y, z=d$z, 
        type="scatter3d",
        mode="markers",
        color=d$isSeed,
        colors = c("#00FF00", "#FF0000")) %>% layout(showlegend = FALSE)


# Check dispersal works for all combinations of cells and proportions
combs = expand.grid("props" = seq(0, 1, 0.05),
                    "cells" = seq(2, 200, 1))

# Test that all combinations can complete without stalling in dispersal finding
invisible(mapply(create.blastocyst, n.cells=combs$cells, prop.aneuploid=combs$props, dispersion=1))
