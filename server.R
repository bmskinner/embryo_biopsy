# Control the server
library(shiny)
library(shinythemes)
library(plotly)
source("fibonacci.R")

# Define the server
function(input, output, session){
  
  
  calculateData = reactive({
    create.blastocyst(n.cells=input$n.cells, 
                          prop.aneuploid = input$proportion, 
                          dispersion = input$dispersal)
  })
  
  output$biopsyPlot = renderPlotly({
    d = calculateData()
    plot_ly(x=d$x, y=d$y, z=d$z, 
            type="scatter3d",
            mode="markers",
            color=d$isSeed,
            colors = c("#00FF00", "#FF0000")) %>% layout(showlegend = FALSE)
  })
  
  output$iterationSummary = renderPlot({

    d = calculateData()
    result = make.samples(d, input$n.samples)
    
    n.euploids = length(result[result==0])
    n.aneuploids  = length(result) - n.euploids
    ratio = (n.euploids / length(result))*100
    
    # Plot the results
    hist(result, xlim=c(-0.5,input$n.samples+0.5),
         breaks=seq(-0.5, input$n.samples+0.5, 1),
         xlab = "Number of aneuploid cells in biopsy",
         main = paste("Biopsying",input$n.samples, 
                      "cells from this blastocyst \nwould give only euploid cells in", 
                      format(ratio, nsmall=1, digits = 3), "% of biopsies"))
  })


}
