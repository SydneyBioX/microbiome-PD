library(shiny)
library(phyloseq)
library(plotly)
library(propr)
library(RColorBrewer)
source("./utils.R")

load(file = "./corr.RData")
ps <- data$ps
clinical.df <- data$clinical.df
vartype <- data$vartype
PD <- as.factor(clinical.df[,1])

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
#    observeEvent(input$plot, {
#            Microbiome <<- microtrans(ps, input$trans, input$Levels, input$data)
##            Clinical <<- clinicaltrans(clinical.df, vartype, input$data)
#            R <<- cal_corr(Microbiome, Clinical, vartype, input$D1, input$D2, input$trans, input$corr.method, input$data)
#            OTU_name <<- OTU_info(ps)
#        })
    
    observeEvent(input$plot,{
        output$heatmap <- renderPlotly({
            Microbiome <<- microtrans(ps, input$trans, input$Levels, input$data, PD)
            Clinical <<- clinicaltrans(clinical.df, vartype, input$data, PD)
            R <<- cal_corr(Microbiome, Clinical, vartype, input$D1, input$D2, input$trans, input$corr.method, input$data, PD)
            OTU_name <<- OTU_info(ps)
            p <- try(plot_heat(R, input$r0, input$mode))
            if('try-error' %in% class(p))
            {
                ggplotly()
            }else{
                p
            }
        })
    })
    
    output$scatterplot <- renderPlotly({
        event.data <- event_data(event = "plotly_click", source = "heatmap")
        if(!is.null(event.data)){
            p <- try(plot_scatter(Microbiome, Clinical, R, input$r0, input$mode, event.data$x, event.data$y, input$D1, input$D2, input$data, PD), silent = T)
            if('try-error' %in% class(p))
            {
                ggplotly()
            }else{
                p
            }
            
        }
        
    })
    output$otuinfo <- renderDataTable({
        OTUinfo <- OTU_info(ps)
        DT::datatable(OTUinfo)
    })
    
    output$orderedby <- renderUI({
        ps_df <<- RelAbundChart(ps, input$Levels1)
        selectInput("ordered","ordered by:",choices = unique(levels(ps_df[,1])))
    })
    
    observeEvent(input$plot1,{
        output$sampleplot <- renderPlotly({
            p <- sample_plot(ps_df,"Sample plot",input$ordered, input$Levels1)
        })
    }
        
    )
   

})
