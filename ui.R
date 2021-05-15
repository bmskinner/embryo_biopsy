# Control the UI
library(shiny)
library(shinythemes)

# Define the UI
fluidPage(theme = shinytheme("lumen"),
            titlePanel("Aneuploidies seen in biopsies"),
            sidebarLayout(
              sidebarPanel(
                numericInput(inputId = "dispersal",
                             label = strong("Dispersal of aneuploid cells (0-1)"),
                             value = 0.1,
                             min=0,
                             max = 1,
                             step = 0.05),
                
                numericInput(inputId = "proportion",
                             label = strong("Proportion of aneuploid cells (0-1)"),
                             value = 0.1,
                             min=0,
                             max = 1,
                             step = 0.05),
                
                numericInput(inputId = "iterations",
                             label = strong("Simulations to run (1-100)"),
                             value = 1,
                             min=1,
                             max = 100,
                             step = 1),
                
                numericInput(inputId = "dim.x",
                             label = strong("Number of columns"),
                             value = 11,
                             min=10,
                             max = 100,
                             step = 1),
                
                numericInput(inputId = "dim.y",
                             label = strong("Number of rows"),
                             value = 11,
                             min=10,
                             max = 100,
                             step = 1)
               
                
              ),
              mainPanel(
                plotOutput(outputId = "biopsyPlot", height = "300px"),
                plotOutput(outputId = "iterationSummary", height = "300px")
              ))
)