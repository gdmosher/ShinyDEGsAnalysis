## app.R ##
library(shiny)
library(shinydashboard)
library(systemPipeR)

# library(systemPipeRdata)
# library(tidyverse)
# library(stringr) # for str_locate() in arendt shell-app
# # install.packages("rhandsontable")
# library(rhandsontable)
# # install.packages("ssh")
# library(ssh)
# # install.packages("shinyAce")
# library(shinyAce)
# # install.packages("shinyjs")
# library(shinyjs)
# # install.packages("V8")
# library(V8)
# library(DT)

source("module-DEGsAnalysis.R")
workDir <<- "sda" # global for use in moduleDEGsAnalysis ui

ui <- dashboardPage(
  dashboardHeader(title = "ShinyDEGsAnalysis"),
  ## Sidebar content
  dashboardSidebar(
    sidebarMenu(
      menuItem("Analysis of DEGs", tabName = "DEGsAnalysisInput6-module", icon = icon("dashboard")),
      menuItem("Clustering and heat maps", tabName = "DEGsAnalysisInput8-module", icon = icon("dashboard"))
    )
  ),
  ## Body content
  dashboardBody(
    tabItems(
      tabItem(tabName = "DEGsAnalysisInput6-module",
              module.DEGsAnalysisInput6(id="module.DEGsAnalysis_2ndInstance"), p1=''), # id here must match callModule(id) below
      tabItem(tabName = "DEGsAnalysisInput8-module",
              module.DEGsAnalysisInput8(id="module.DEGsAnalysis_2ndInstance"), p1='') # id here must match callModule(id) below
      
    )
  )
)

server <- function(input, output) {

  module.DEGsAnalysis <- callModule(module.DEGsAnalysis, "module.DEGsAnalysis_2ndInstance", "")

  # set.seed(122)
  # histdata <- rnorm(500)
  # 
  # output$plot1 <- renderPlot({
  #   data <- histdata[seq_len(input$slider)]
  #   hist(data)
  # })
}

shinyApp(ui, server)