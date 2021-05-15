# Shiny app
library(tidyverse)
library(shiny)
library(shinythemes)

# The modelling is handled in a separate file
source("grid.R")
source("ui.R")
source("server.R")

shinyApp(ui = create.ui(), server = server)
