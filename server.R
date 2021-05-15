# Control the server
library(shiny)
library(shinythemes)
source("grid.R")

# Define the server
function(input, output, session){
  calculateData = reactive({

    result = lapply(1:input$iterations, make.count, 
                              prop.aneuploids= input$proportion, 
                              dispersion = input$dispersal, 
                              make.chart=T,
                              dim.x = input$dim.x,
                              dim.y = input$dim.y)
    
    
    
  })
  
  output$biopsyPlot = renderPlot({
    calculateData()[[1]][['plot']]
  })
  
  output$iterationSummary = renderPlot({
    # Extract counts
    counts.and.plots = calculateData()
    counts = sapply(1:input$iterations, function(i) counts.and.plots[[i]][['n.euploid']])
    # Plot the results
    hist(counts, xlim=c(-0.5,5.5), 
         breaks=c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5),
         xlab = "Number of euploid cells in biopsy",
         main = paste("Results from", input$iterations, "simulations"))
  })


}
