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
source("shiny_hackathon/R/plotGfa.R")
source("shiny_hackathon/R/plot_arrow_linear_annotations.R")

# functions
# function to change rgb input into hex color
rgb_to_hex <- function(rgb_comm){
  rgb_comm = strsplit(split = ",", x = rgb_comm) %>% unlist()
  return(rgb(red = rgb_comm[1], green = rgb_comm[2], blue = rgb_comm[3], maxColorValue = 255))
}


# User interface ----
ui <- fluidPage(
   tags$head(
      tags$style("html, body { height: 100%; width: 100%}"),
      tags$style("#panel1 {height: 100px; position: fixed}"),
      tags$style("#panel2 {
              overflow: auto;
              background: orange;
              margin-left: width:20%;
          }")
      ),
		 ## Title 
		titlePanel("GFA Visualization"),
				## Sidebar content
				sidebarPanel(
						# Input: Select a file ----
						## Input rGFA file
						fileInput(
							"rgfaFile","rGFA data input",
							multiple=T, 
							buttonLabel="Browse",
							placeholder="No file selected",
							accept = c("text/plain",".txt")),

						## Input BED file
						fileInput(
							"bed",
						  	"BED Annotation",
							multiple=F, 
							buttonLabel="browse",
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
					),
       absolutePanel(id = "panel2", 
                     top = "50%", left = "35%", height = "40%", width = "60%", right = "10%",bottom = "10%",
                     fluidRow(## Plot Ouput of linear visualization
                       p("Linear Visualization Window"),
                       uiOutput("bed_plots.ui"),
                     ),
                 
    ),
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

## Linear Annotations
	  file_list = reactiveValues()
    observe({
       if( !is.null(input$bed) ){ # & !(input$bed %in% file_list$dList) ){
         file_list$dList = append( isolate(file_list$dList) , isolate(input$bed$datapath) )
       }
     })
    output$file_list <- renderPrint({
      req(input$bed)
      paste(file_list$dList, sep = ",")
     })

    observe({
      file_list
    })
    ## Linear visualization
    #p = NULL
    get_bed_df <- reactive({
      req(input$bed)
      bed_df = load_annotation_bed(bed_path = input$bed$datapath)
      bed_df
    })
    plotHeight <- reactive({
      cur_bed = get_bed_df()
      100 * length( unique(cur_bed$contig) )
    }
    ) 
    output$bed_plots <- renderPlot( {
      req(get_bed_df())
      if (is.null(get_bed_df()) ){
        plot.new()
        return()
      }
      cur_df = load_annotation_bed(bed_path = input$bed$datapath, color_col = 9)
      p = plot_bed_annot_track(track_name = "testing", bed_df = cur_df, p = NULL)
      p
    })
    output$bed_plots.ui <- renderUI({
      plotOutput("bed_plots", 
                 height = plotHeight(), 
                 click = "plot_click",
                 dblclick = "plot_dblclick",
                 hover = "plot_hover",
                 brush = "plot_brush",
                 inline = F)
    })

	## Linear brushing visualization
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
		res <- brushedPoints(bed_df(),input$plot_brush,xvar = "start",yvar = NULL,allRows = FALSE)
		if (nrow(res) == 0)
			return()
		res[c("contig","start","stop","name","strand")]
	})
}

# Run app ----
shinyApp(ui, server)
