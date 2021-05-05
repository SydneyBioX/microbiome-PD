pkgs<- c("shiny","shinyjs","shinythemes","ggplot2","plotly","phyloseq","reshape2","tidyverse",
         "DT","propr")
for (i in 1:length(pkgs)){
    pkg <- pkgs[i]
    test <- require(pkg,character.only = T)
    if(!test){
        install.packages(pkg,ask = F,update = F)
        require(pkg,character.only = T)
    }
}


library(shiny)

shinyUI(fluidPage(
    navbarPage(
        "Cross-Sectional Analysis",
        tabPanel("Correlation Analysis",
                 sidebarLayout(
                     sidebarPanel(
                         selectInput("D1",
                                     "Input data1:",
                                     choices = c("Microbiome","Clinical"),
                                     selected = "Microbiome"),
                         selectInput("D2",
                                     "Input data2:",
                                     choices = c("Microbiome","Clinical"),
                                     selected = "Microbiome"),
                         selectInput("Levels",
                                     "Taxa Levels for microbiome data:",
                                     choices = c("Phylum","Order","Family","Genus","ASV"),
                                     selected = "Genus"),
                         sliderInput("r0",
                                     label = "Correlation threshold:",
                                     min = 0,
                                     max = 1,
                                     value = 0.5,
                                     step = 0.05),
                         radioButtons("trans","Microbiome data transformation:",choices = c("Count","Composition","Log"), selected = "Composition"),
                         radioButtons("mode","Output values:",choices = c("coefficient","-1-0-1"),selected = "coefficient"),
                         radioButtons("corr.method","Correlation method:",choices = c("pearson","spearman"), selected = "spearman"),
                         radioButtons("data","data included:", choices = c("All", "PD"), selected = "All"),
                         actionButton("plot", "Analysis")
                         
                     ),
                     
                     mainPanel(
                         plotlyOutput("heatmap"),
                         plotlyOutput("scatterplot")
                     ))),
        
        tabPanel("ASV info",
                 mainPanel(DT::dataTableOutput("otuinfo"))),
        
        tabPanel("Sample plot",
                 sidebarLayout(
                     sidebarPanel(
                         selectInput("Levels1",
                                     "Taxa Levels for microbiome data:",
                                     choices = c("Phylum","Order","Family","Genus"),
                                     selected = "Genus"),
                         uiOutput("orderedby"),
                         actionButton("plot1","Analysis")
                         
                     ),
                     mainPanel(
                         plotlyOutput("sampleplot")
                     )
                 ))
    )

    
))
