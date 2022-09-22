# Load packages ----
library(shiny)
library(ggplot2)
library(ggforce)
library(GenomicRanges)
library(gggenes)
library(tidyverse)
library(readr)


# Source helper functions -----
# source("helpers.R")
loadSupport()
#source("shiny_hackathon/R/plotGfa.R")
#source("shiny_hackathon/R/plot_arrow_linear_annotations.R")

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
						# Not so sure about this part yet
						selectInput(
						  "select_graph","Graph",selected = NULL, multiple = FALSE,list("outputfile1","outputfile2")),
						
						
						## Download Fasta file and Bed file
						downloadButton("Fasta_download", "FASTA file Download",icon = shiny::icon("download")),
						downloadButton("BED_download", "BED file Download",icon = shiny::icon("download"))
						),
				absolutePanel(id = "panel1",
				             top = "5%", left = "35%", height = "40%", width = "60%", right = "10%",
				             plotOutput(
				               "ggdag",
				               "Graph Visualization Window",
				               width="100%",
				               height="100px",
				               brush=brushOpts(id="plot_brush")
				             )
				             ),
				# mainPanel(
				# 	        "Graph Visualization Window",
				# 
				# 	        ## Plot Ouput of graphic visualization
				# 	        plotOutput(
				# 			"ggdag",
				# 			"Graph Visualization Window",
				# 			width="100%",
				# 			height="100px",
				# 			brush=brushOpts(id="plot_brush")
				# 		        ),
				# 	),
       absolutePanel(id = "panel2", 
                     top = "50%", left = "35%", height = "30%", width = "60%", right = "5%",bottom = "10%",
                     fluidRow(## Plot Ouput of linear visualization
                       p("Linear Visualization Window"),
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
server <- function(input, output) {
  spacer.width <- reactive({
    input$spacer
    })
	# Read the rGFA file
	output$contents <- renderTable({
		req(input$rgfaFile)
	  df <- readGfa(gfa.file = input$rgfaFile$datapath, 
	                store.fasta = 'FALSE')
		return(df$segments)
	})

	# Plot the graph
	output$ggdag <- renderPlot({
	  req(input$rgfaFile)
	  df <- readGfa(gfa.file = input$rgfaFile$datapath, 
	                store.fasta = 'FALSE')
	  segms.gr <- GRanges(seqnames = 'nodes', ranges = IRanges(start = 1, end = df$segments$LN), id= df$segments$segment.id)
	  shifts <- width(segms.gr)
	  shifts <- cumsum(shifts + spacer.width())
	  segms.gr[-1] <- shift(segms.gr[-1], shift = shifts[-length(shifts)])
	  segms.df <- as.data.frame(segms.gr)
	  nodes <- data.frame(x=c(rbind(segms.df$start, segms.df$end)),
	                      y=0,
	                      group=rep(1:nrow(segms.df), each=2))
	 
	  ## Plt nodes/segments
	  node.plt <- ggplot(nodes, aes(x = x, y = y, group=group)) +
	    geom_shape(radius = unit(0.5, 'cm'))
	  
	  ## Define links
	  links <- data.frame(from=as.numeric(gsub(df$links$from, pattern = 's', replacement = '')),
	                      to=as.numeric(gsub(df$links$to, pattern = 's', replacement = '')))
	  arcs <- c(rbind(segms.df$end[links$from], segms.df$start[links$to]))
	  
	  ## Make geom_curve
	  arcs.df <- data.frame(x=arcs[c(TRUE, FALSE)],
	                        xend=arcs[c(FALSE, TRUE)],
	                        y=0,
	                        yend=0,
	                        shape=ifelse((links$to - links$from) == 1, 'line', 'curve'))
	  
	  curve.graph.plt <- node.plt + 
	    geom_segment(data=arcs.df[arcs.df$shape == 'line',], aes(x = x, y = y, xend=xend, yend=yend), arrow = arrow(length = unit(0.01, "npc")), inherit.aes = FALSE) +
	    geom_curve(data=arcs.df[arcs.df$shape == 'curve',], aes(x = x, y = y, xend=xend, yend=yend), ncp = 100, curvature = -1, arrow = arrow(length = unit(0.01, "npc")), inherit.aes = FALSE) +
	    geom_point(data=arcs.df,  aes(x = x, y = y), inherit.aes = FALSE)
	  
	  ## Make geom_bezier
	  levels <- (links$to - links$from) - 1
	  arcs.df <- data.frame(x=rep(arcs, each=2),
	                        y=do.call(c, lapply(levels, function(x) c(0,x,x,0))),
	                        group=rep(1:nrow(links), each=4))
	  
	  node.plt + 
	    geom_bezier(data=arcs.df, aes(x = x, y = y, group=group), arrow = arrow(length = unit(0.01, "npc")), inherit.aes = FALSE)
	  })
					
	  output$info <- renderText({
		  xy_str <- function(e) {
			  if(is.null(e)) return("NULL\n")
			  paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
		  }
		  xy_range_str <- function(e) {
			  if(is.null(e)) return("NULL\n")
			  paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1), 
				 " ymin=", round(e$ymin, 1), " ymax=", round(e$ymax, 1))
		  }
		  paste0(
			  "click: ", xy_str(input$plot_click),
			  "dblclick: ", xy_str(input$plot_dblclick),
			  "hover: ", xy_str(input$plot_hover),
			  "brush: ", xy_range_str(input$plot_brush)
		  )
	  })
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
}

# Run app ----
shinyApp(ui, server)
