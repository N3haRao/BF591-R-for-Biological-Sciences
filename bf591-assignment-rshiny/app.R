## Author: Neha Rao
## neharao@bu.edu
## BU BF591
## Assignment 7

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker) # you might need to install this

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Neha Rao: BF591 - Assignment 7"),
  sidebarLayout(
    sidebarPanel(
      fileInput("exp_file", label = "Load differential expression results:", accept = c(".csv")),
      
      p('A volcano plot can be generated with "log2 fold-change" on the x-axis and "p-adjusted" on the y-axis.'),
      
      radioButtons("x_axis", "Choose the column for the X-axis:", 'log2FoldChange',
                   choices = c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj')),
      
      radioButtons("y_axis", "Choose the column for the Y-axis:", 'padj',
                   choices = c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj')),
      
      colourInput("base_point_col", 'Base point color', "#56e0ff"),
      
      colourInput("highlight_col", 'Highlight point color', "#f00206"),
      
      sliderInput(inputId = "p_magnitude", min = -300, max = 0,
                  label = "Select the magnitude of the p adjusted coloring:", value = -150, step = 1),
      
      submitButton("Submit", icon("check"))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Table",
                 tableOutput('table')
        ),
        tabPanel("Plot",
                 plotOutput('volcano')
        )
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # Load Data
  load_data <- reactive({
    file <- input$exp_file
    if (is.null(file)) {
      return(NULL)
    }
    return(read.csv(file$datapath, header = TRUE))
  })
  
  # Volcano plot
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    if (is.null(dataf)) {
      return(NULL)
    }
    y_label <- sprintf("-log10(%s)", y_name)
    col_label <- sprintf("padj < (1*(10^%s))", slider)
    
    vol_plot <- ggplot(dataf, aes(x = get(x_name), y = -log10(get(y_name)), color = padj < (1*(10^slider)))) +
      geom_point() +
      scale_colour_manual(name = col_label, values = setNames(c(color1, color2, 'black'), c(T, F, NA))) +
      xlab(x_name) +
      ylab(y_label) +
      theme_linedraw()
    
    return(vol_plot)
  }
  
  # Draw and filter table
  draw_table <- function(dataf, slider) {
    if (is.null(dataf)) {
      return(NULL)
    }
    
    dataf_filtered <- dplyr::filter(dataf, padj < (1 * (10^slider)))
    dataf_filtered$pvalue <- lapply(dataf_filtered$pvalue, formatC, digits = 9)
    dataf_filtered$padj <- lapply(dataf_filtered$padj, formatC, digits = 9)
    
    return(dataf_filtered)
  }
  
  output$volcano <- renderPlot({
    dataf <- load_data()
    x_name <- input$x_axis
    y_name <- input$y_axis
    slider <- input$p_magnitude
    color1 <- input$base_point_col
    color2 <- input$highlight_col
    
    p <- volcano_plot(dataf, x_name, y_name, slider, color1, color2)
    return(p)
  })
  
  output$table <- renderTable({
    dataf <- load_data()
    slider <- input$p_magnitude
    
    dataf_filtered <- draw_table(dataf, slider)
    return(dataf_filtered)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
