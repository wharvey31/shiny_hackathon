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
							multiple=T, 
							buttonLabel="Browse",
							placeholder="No file selected",
							accept = c("text/plain",".txt")),
						## Input BED file
						fileInput(
							"bed",
						  	"BED file input",
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
						  "select_graph","Graph",selected = NULL, multiple = FALSE,list("outputfile1","outputfile2")),
						
						
						## Download Fasta file and Bed file
						downloadButton("Fasta_download", "FASTA file Download",icon = shiny::icon("download")),
						downloadButton("BED_download", "BED file Download",icon = shiny::icon("download"))
						),
				
				mainPanel(
					        "Graph Visualization Window",

					        ## Plot Ouput of graphic visualization
					        plotOutput(
							"ggdag",
							"Graph Visualization Window",
							width="100%",
							height="400px",
							brush=brushOpts(id="plot_brush")
						        ),

					        ## Plot Ouput of linear visualization
					        "Linear Visualization Window",
					        plotOutput(
							"bed_plots",width="100%",height="100px",click = "plot_click",
							dblclick = "plot_dblclick",
							hover = "plot_hover",
							brush = "plot_brush"
						        ),

					        ## Ouput of user graphical interactive information
					        verbatimTextOutput("info"),
					        tableOutput("contents"),
							"File List",
							verbatimTextOutput('file_list')
					)
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
					  
	  ## Linear visualization
	  output$bed_plots <- renderPlot({
		  req(input$bed)
		  bed_df <- read_tsv(input$bed$datapath, col_names = T )
		  colnames(bed_df) <- c("contig","start","stop","name",".","strand","..","...","rgb")
		  rgb_to_hex <- function(rgb_comm){
			  rgb_comm = strsplit(split = ",", x = rgb_comm) %>% unlist()
			  return(rgb(red = rgb_comm[1], green = rgb_comm[2], blue = rgb_comm[3], maxColorValue = 255))
		  	  }
		  bed_df$hex_color<-sapply(bed_df$rgb, FUN = rgb_to_hex)

		  ggplot(data = bed_df) +
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
}

# Run app ----
shinyApp(ui, server)
