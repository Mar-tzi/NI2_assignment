#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
# app.R
library(shiny)
library(dplyr)
library(ggplot2)
library(colourpicker)  # For color picker
library(DT)  # For rendering tables

ui <- fluidPage(
  titlePanel("Gene Expression Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV File", accept = c(".csv")),
      selectInput("gene", "Select Gene to Analyze:", choices = NULL),
      selectInput("controlSample", "Select Control (Basal) Sample:", choices = NULL),
      selectInput("removeSample", "Select Sample to Remove (Optional):", choices = c("None")),
      colourInput("plotColor", "Choose Plot Color", value = "#1f77b4"),  # Default color
      textInput("plotTitle", "Plot Title", value = "Gene Expression Plot"),
      textInput("xLabel", "X Axis Label", value = "Sample"),
      textInput("yLabel", "Y Axis Label", value = "Fold Change"),
      selectInput("plotTheme", "Select Plot Theme:", choices = c("Minimal" = "minimal", "Classic" = "classic", "Dark" = "dark")),
      sliderInput("barWidth", "Bar Width:", min = 0.1, max = 1, value = 0.7, step = 0.1)
    ),
    
    mainPanel(
      plotOutput("genePlot", height = "600px"),
      DTOutput("sampleTable")
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive expression to read and process data
  data <- reactive({
    req(input$file)
    df <- read.csv(input$file$datapath, stringsAsFactors = FALSE)
    
    # Calculate average Ct values per replicate per sample
    avg_df <- df %>%
      group_by(Gene, Sample) %>%
      summarise(Avg_Ct = mean(Ct), .groups = 'drop')
    
    # Calculate average Ct for beta-actin per sample
    beta_actin_avg <- avg_df %>%
      filter(Gene == "beta_actin") %>%
      group_by(Sample) %>%
      summarise(Beta_Actin_Avg = mean(Avg_Ct), .groups = 'drop')
    
    # Normalize Ct values by beta-actin
    normalized_df <- avg_df %>%
      filter(Gene != "beta_actin") %>%
      left_join(beta_actin_avg, by = "Sample") %>%
      mutate(Delta_Ct = Avg_Ct - Beta_Actin_Avg)
    
    # Calculate Delta-Delta Ct and Fold Change
    delta_delta_df <- normalized_df %>%
      group_by(Gene) %>%
      mutate(
        Mean_Delta_Ct_Control = mean(Delta_Ct[Sample == input$controlSample], na.rm = TRUE),
        Delta_Delta_Ct = Delta_Ct - Mean_Delta_Ct_Control,
        Fold_Change = 2^(-Delta_Delta_Ct)
      ) %>%
      ungroup()
    
    # Ensure the gene and control samples are updated in the input widgets
    updateSelectInput(session, "gene", choices = unique(delta_delta_df$Gene))
    updateSelectInput(session, "controlSample", choices = unique(df$Sample))
    updateSelectInput(session, "removeSample", choices = c("None", unique(df$Sample)))
    
    delta_delta_df
  })
  
  output$genePlot <- renderPlot({
    df <- data()
    if (is.null(df)) return(NULL)
    
    selected_gene <- input$gene
    plot_color <- input$plotColor  # Get selected plot color
    plot_title <- input$plotTitle
    x_label <- input$xLabel
    y_label <- input$yLabel
    plot_theme <- switch(input$plotTheme,
                         "minimal" = theme_minimal(),
                         "classic" = theme_classic(),
                         "dark" = theme_dark())
    bar_width <- input$barWidth
    
    filtered_data <- df %>%
      filter(Gene == selected_gene)
    
    # Apply the sample removal if specified
    if (input$removeSample != "None") {
      filtered_data <- filtered_data %>%
        filter(Sample != input$removeSample)
    }
    
    if (nrow(filtered_data) == 0) {
      showNotification("No data available for the selected gene.", type = "error")
      return(NULL)
    }
    
    ggplot(filtered_data, aes(x = Sample, y = Fold_Change, fill = Sample)) +
      scale_fill_manual(values = rep(plot_color, length(unique(filtered_data$Sample)))) +
      labs(title = plot_title,
           x = x_label, y = y_label) +
      plot_theme +
      geom_bar(stat = "identity", position = position_dodge(), width = bar_width)
  })
  
  output$sampleTable <- renderDT({
    df <- data()
    if (is.null(df)) return(NULL)
    
    # Generate the table with default settings, summarizing by sample
    df %>%
      filter(Gene == input$gene) %>%
      filter(Sample != input$removeSample) %>%
      group_by(Sample) %>%
      summarise(Fold_Change = mean(Fold_Change), .groups = 'drop') %>%
      datatable(options = list(pageLength = 10, dom = 't'))
  })
}

shinyApp(ui = ui, server = server)



