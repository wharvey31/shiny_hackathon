# Load packages ----
library(shiny)
library(ggplot2)
library(ggforce)
library(GenomicRanges)
library(gggenes)
library(tidyverse)
library(readr)
# library(Cairo)


# Source helper functions -----
# source("helpers.R")
loadSupport()

# functions
# function to change rgb input into hex color
rgb_to_hex <- function(rgb_comm){
  rgb_comm = strsplit(split = ",", x = rgb_comm) %>% unlist()
  return(rgb(red = rgb_comm[1], green = rgb_comm[2], blue = rgb_comm[3], maxColorValue = 255))
}


# User interface ----
ui <- fluidPage(
		 ## Title 
		titlePanel("GFA Visualization"),
				## Sidebar content
				sidebarPanel(
						# Input: Select a file ----
						## Input rGFA file
						fileInput(
							"rgfaFile","rGFA data input",
							multiple=FALSE, 
							buttonLabel="Browse",
							placeholder="No file selected",
							accept = c("text/plain",".txt")),
						 ## Input BED file
						 fileInput(
							 "bed",
							 "BED file input",
							 multiple=FALSE,
							 buttonLabel="Browse",
							 placeholder="No file selected",
							 accept = c("text/plain",".txt")),
						## Input GAF file
						fileInput(
							"GAF_input2",
							"GAF data input",
							multiple=FALSE, 
							buttonLabel="Browse",
							placeholder="No file selected"
							),
						## Input Annotation file
						fileInput(
							"Annotation_Input",
							"Add Annotation",
							multiple=FALSE, 
							buttonLabel="Browse",
							placeholder="No file selected"
							),
						# Not so sure about this part yet
						selectInput(
						  "select_graph",
						  "Graph",
						  selected = NULL, 
						  multiple = FALSE,
						  list("outputfile1","outputfile2")
						  ),
						## Download Fasta file and Bed file
						downloadButton("Fasta_download", 
						               "FASTA file Download",
						               icon = shiny::icon("download")
						               ),
						downloadButton("BED_download", 
						               "BED file Download",
						               icon = shiny::icon("download"))
						),
				
				mainPanel(
					        # "Graph Visualization Window",

					        ## Plot Ouput of graphic visualization
					        plotOutput(
      							"ggdag",
      							"Graph Visualization Window",
      							width="100%",
      							height="400px",
      							brush=brushOpts(id="graph_brush", resetOnNew = TRUE),
							      click="graph_click",
      							dblclick = "graph_dblclick",
						        ),
					        tableOutput("graph_point"),
					        ## Plot Ouput of linear visualization
					        # Linear Visualization Window
					        h3("Linear Visualization Window"),
						plotOutput(
							"bed_plots",width="100%",height="100px",click = "plot_click",
							dblclick = "plot_dblclick",
							hover = "plot_hover",
							brush = brushOpts(id = "plot_brush",direction="x")
							),
					        ## Ouput of user graphical interactive information
						h3("Brushed info"),
						verbatimTextOutput("info"),
						h4("Brushed gene info"),
						tableOutput("plot_brushedpoints"),
					)
		)

# Server logic ----
server <- function(input, output) {
	# Read the rGFA file
	ranges <- reactiveValues(x = NULL, y = NULL)
	observeEvent(input$graph_dblclick, {
		brush <- input$graph_brush
		if (!is.null(brush)) {
			ranges$x <- c(brush$xmin, brush$xmax)
		} else {
			ranges$x <- NULL
		}
	})
	graph_df <- reactive({
		readGfa(gfa.file = input$rgfaFile$datapath, 
			store.fasta = 'FALSE')
	})
	# Plot the graph
	output$ggdag <- renderPlot({
		# Wait for rgfa input
		req(input$rgfaFile)
		# Plot df and add cartesian 
		plotGfa(gfa.tbl=graph_df()) + coord_cartesian(xlim = ranges$x, ylim = NULL, expand = FALSE)
	})
	output$graph_point <- renderTable({
		req(input$graph_click)
		brushedPoints(graph_df()$segments, input$graph_brush, xvar = "SO", yvar = "SR")
	})
	## Linear visualization
	bed_df <- reactive({
		data_tsv <- read_tsv(input$bed$datapath, col_names = T )
		colnames(data_tsv) <- c("contig","start","stop","name",".","strand","..","...","rgb")
		data_tsv$hex_color <- sapply(data_tsv$rgb, FUN = rgb_to_hex)
		data_tsv
	})
	
	output$bed_plots <- renderPlot({
		req(input$bed)
		ggplot(data = bed_df()) +
		gggenes::geom_gene_arrow(mapping =  aes(xmin = start, xmax = stop, y = 1, fill = hex_color)) +
		theme(legend.position="none")
	})
  	
	output$info <- renderText({
		xy_str <- function(e) {
			if(is.null(e)) return("NULL\n")
			paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
		}
		
		xy_range_str <- function(e) {
			if(is.null(e)) return("NULL\n")
			paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1))
		}
		
		paste0(
			"click: ", xy_str(input$plot_click),
			"dblclick: ", xy_str(input$plot_dblclick),
			"hover: ", xy_str(input$plot_hover),
			"brush: ", xy_range_str(input$plot_brush)
		)
	})
	output$plot_brushedpoints <- renderTable({
		bed_df <- read_tsv(input$bed$datapath, col_names = T )
		colnames(bed_df) <- c("contig","start","stop","name",".","strand","..","...","rgb")
		res <- brushedPoints(bed_df,input$plot_brush,xvar = "start",yvar = NULL,allRows = FALSE)
		if (nrow(res) == 0)
			return()
		res[c("contig","start","stop","name","strand")]
	})
}

# Run app ----
shinyApp(ui, server)
