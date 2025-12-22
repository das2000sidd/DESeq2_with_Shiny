#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(DESeq2)
library(DT)
library(ggplot2)
library(org.Hs.eg.db)
library(ggrepel)
library(grid)
library(dplyr)


# Define server logic required to draw a histogram

ui <- fluidPage(
  titlePanel("Bulk RNA-seq differential expression visualisation for DESeq2 Results"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("de_file", "Upload DESeq2 results (.csv or .txt)",
                accept = c(".csv", ".txt", ".tsv")),
      
      radioButtons(
        "sep",
        "Separator",
        choices = c(Comma = ",", Tab = "\t"),
        selected = "\t"
      ),
      
      radioButtons(
        "org",
        "Species",
        choices = c(Human = "Homo sapiens", Mouse = "Mus musculus"),
        selected = "Homo sapiens"
      ),
      
      tags$hr(),
      sliderInput("padj", "Adjusted p-value cutoff:",
                  value = 0.05, min = 0, max = 1, step = 0.01),
      sliderInput("log2FoldChange", "Log2 fold change cutoff:",
                  value = 1, min = 0, max = 100, step = 0.1)
    ),
    
    mainPanel(
      tabsetPanel(
        #tabPanel("Summary Table", DTOutput("table")),
        #tabPanel("Volcano Plot", plotOutput("volcano")),
        #tabPanel("MA Plot", plotOutput("maplot"))
        ## Create a tab panel
        ## DTOutput - create a container for table
        ## plotOutput - Create an plot or image output element
        ## tabPanel - If we want to display results, it creates a panel for each one of them
        tabPanel("Results Table",
                 DTOutput("results_table"),
                 br(),
                 downloadButton("download_table", "Download DE Table (CSV)")
        ),
        
        tabPanel("Volcano",
                 plotOutput("volcano"),
                 br(),
                 downloadButton("download_volcano", "Download Volcano Plot")
        ),
        
        tabPanel("MA plot",
                plotOutput("maplot"),
                 br(),
                 downloadButton("download_maplot", "Download MA Plot")
        ),
        
      )
    )
  )
)  

