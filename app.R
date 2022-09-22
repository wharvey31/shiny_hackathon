# Load packages ----
library(shiny)
library(ggplot2)
library(ggforce)
library(GenomicRanges)
library(gggenes)
library(tidyverse)
library(readr)
library(shinyjs)
library(wesanderson)
# library(Cairo)


# Source helper functions -----

loadSupport()

# functions
# function to change rgb input into hex color

# User interface ----
ui <- fluidPage(
	  tags$head(
      tags$style("html, body { height: 100%; width: 100%}"),
      tags$style("#title_panel {
                 background : blue;
                 margin-left: width:20%"
      ),
      tags$style("#side_panel {
                 background : green;
                 margin-left: width:20%"
                 ),
      tags$style("#panel1 {
      background: red;
                 margin-left: width:100%;
                 }"),
      tags$style("#panel2 {
              overflow: auto;
              background: orange;
              margin-left: width:100%;
          }"),
      tags$style("#panel3 {
              overflow: auto;
              background: purple;
              margin-left: width:20%;
          }"),
      ),
		 ## Title 
		absolutePanel(id = "title_panel",
		              top = "0%", left = "1%", height = "15%", width = "30%", bottom = "90%",
		              h3("GFA Visualization") ),
				## Sidebar content
				absolutePanel(id = "side_panel",
				              top = "20%", left = "1%", height = "70%", width = "30%", bottom = "20%", #right = "80%"
					    # Initialize shiny jus
						shinyjs::useShinyjs(),
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
						uiOutput("contig_selection"),
						# Haplotype Selection
						shinyjs::disabled(
							selectInput(
								"select_graph",
								"Graph",
								selected = NULL, 
								multiple = FALSE,
								list()
							)
						),
						shinyjs::disabled(
							radioButtons(
								"showFrequency",
								"Show Frequency:",
								c("None" = "none",
									"Width" = "width",
									"Color" = "color"),
								inline=TRUE
							)
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

				absolutePanel(id = "panel1",
				             top = "5%", left = "35%", height = "40%", width = "60%", right = "10%",
				             plotOutput(
				               "ggdag",
				               "Graph Visualization Window",
				               width="100%",
				               height="100%",
				               brush=brushOpts(id="graph_brush", resetOnNew = TRUE),
				               click="graph_click",
				               dblclick="graph_dblclick",
				             )
				             ),
       absolutePanel(id = "panel2", 
                     top = "50%", left = "35%", height = "100%", width = "60%", right = "5%",bottom = "10%",
                     fluidRow(## Plot Ouput of linear visualization
                       # p("Linear Visualization Window"),
                       uiOutput("bed_plots.ui"),
                     ),
                 
    ),
   absolutePanel(id = "panel3", 
                 top = "80%", left = "35%", height = "15%", width = "60%", right = "10%",bottom = "10%",
                 fluidRow(## Plot Ouput of linear visualization
                   p("GAF segments"),
                   uiOutput("GAF_select.ui"),
                 ),
   ),
)

# Server logic ----
server <- function(input, output, session) {
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
      p = plot_bed_annot_track(track_name = "", bed_df = cur_df, p = NULL)
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
    contigs <- reactive({
      cur_bed = get_bed_df()
      get_haplotype_names(cur_bed)
    }
    )
    
    #choices of contig:
    output$contig_selection = renderUI({
      selectInput(inputId = 'contig_selection',label = 'Contig Highlight', contigs())
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

	# Haplotype select population
	haplotypes_df <- reactive({
		gafToLinks(gaf.file = input$GAF_input2$datapath)
	})
	
	haplotypes <- eventReactive(input$GAF_input2,{
		haplotypes_df()$haplotype[!duplicated(haplotypes_df()$haplotype)]
	})
	
	observeEvent(haplotypes(), {
		updateSelectInput(session = session, inputId = "select_graph", choices = haplotypes())
		shinyjs::enable("select_graph")
		shinyjs::enable("showFrequency")
	})

	# Haplotype selection
	haplotype <- eventReactive(input$select_graph, {
		input$select_graph			
	})

	# Show Frequency selection
	show_frequency <- eventReactive(input$showFrequency, {
		showFrequencySelection <- input$showFrequency
		if (showFrequencySelection == "none")
			showFrequencySelection = NULL
		return(showFrequencySelection)
	})

	observeEvent(haplotypes(), {
		output$ggdag  <- renderPlot({
			
			full_plot <- plotGfa(gfa.tbl=graph_df())
			max_abs_value <- max_absolute_value(full_plot$plot_env$arc.height)
			
			haplotype_links <- subset(haplotypes_df(), haplotype==haplotype())
			haplotype_segments <- segments_for_haplotype_links(graph_df()$segments, haplotype_links)
			
			haplotype_info = list(segments = haplotype_segments, links = haplotype_links)
			
			link.frequency <- show_frequency()
			gaf.links <- NULL
			if(!is.null(link.frequency))
				gaf.links <- haplotype_links
			plotGfa(gfa.tbl=haplotype_info, y.limit=max_abs_value, link.frequency=link.frequency, gaf.links=gaf.links) + coord_cartesian(xlim = ranges$x, ylim = NULL, expand = FALSE)
		})
	})
	
}

# Run app ----
shinyApp(ui, server)
