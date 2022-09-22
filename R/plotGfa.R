#' Plot GFA from loaded data table
#' 
#' This function takes an output from readGfa.R function and makes a graph plot (Adjust annotation!!!)
#'
#' @param gfa.tbl A \code{list} of data tables needed for plotting after running 'readGfa' function.
#' @param min.segment.length A minimum size of segment to be plotted.
#' @param min.link.degree A minimum number of each link is traversed in order to be kept.
#' @param spacer.width User defined fraction to the total segment length to be used as node spacer.
#' @param order.by Define a column to be used for node ordering. [TODO]
#' @param layout Overall layout of the graph, either 'linear' or 'circular'.
#' @param shape Either 'rectangle' or 'roundrect' to plot graph nodes.
#' @param arrow.head Default 'closed' TODO what are the other options?
#' @param gaf.links A \code{tibble} table containing links present in each haplotype after alignment to the graph. (GafToLinks output)
#' @param gaf.annotation  A \code{tibble} table containing path coordinates where given annotation map on the graph. (readGaf output)
#' @param link.frequency Visualize link frequency either by the thickness of a link ('width') or a color gradient ('color').
#' @param highlight.haplotype A character string of an haplotype ID to for which the graph path will be highlighted in red color.
#' @importFrom GenomicRanges shift GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors lapply
#' @importFrom ggforce geom_bezier
#' @author David Porubsky, Sean McGee & Karynne Patterson
#' @export
#' 

plotGfa <- function(gfa.tbl=NULL, y.limit=NULL, min.segment.length=0, min.link.degree=0, spacer.width=0.05, order.by='offset', layout='linear', shape='rectangle', arrow.head='closed', gaf.links=NULL, gaf.annotation=NULL, link.frequency=NULL, highlight.haplotype=NULL) {


  ## Check user input ##
  ######################
  ## Get link data from loaded GFA file
  #segments <- gfa.tbl$segments
  links <- gfa.tbl$links
  links.ids <- paste(links$from, links$to, sep = '_')
  
  ## Get link frequency from links if gaf.links are defined
  if (!is.null(gaf.links)) {
    gaf.links.ids <- paste(gaf.links$from, gaf.links$to, sep = '_')
    link.freq <- as.data.frame(table(gaf.links.ids))
    link.freq <- link.freq[order(link.freq$Freq, decreasing = TRUE),]
    links$link.freq <- link.freq$Freq[match(links.ids, link.freq$gaf.links.ids)]
  }
  ## Highlight haplotype
  if (!is.null(highlight.haplotype) & !is.null(gaf.links)) {
    if (highlight.haplotype %in% gaf.links$SN) {
      hap.links <- gaf.links[gaf.links$SN %in% highlight.haplotype,]
      hap.links.ids <- paste(hap.links$from, hap.links$to, sep = '_')
      hap.segments <- unique(hap.links$from, hap.links$to)
    }
  }
  
  ### Prepare data for plotting ###
  #################################
  ## Define segment spacer as fraction of the total length of all segments
  gfa.tbl <- addGraphPlottingCoords(gfa.tbl = gfa.tbl, spacer.width = spacer.width, order.by = order.by)
  segms.df <- gfa.tbl$graph.plt.coords
  
  ## Filter data ##
  #################
  ## Filter segments by size [TODO]
  #segments <- segments[segments$LN >= min.segment.length,]
  #links <- links[links$from %in% segments$segment.id & links$to %in% segments$segment.id,]
  
  ## Filter links by degree
  if (min.link.degree > 0) {
    if (!is.null(gaf.links)) {
      links <- links[links$link.freq > 1,]
      links.segments <- unique(links$from, links$to)
      links.ids <- paste(links$from, links$to, sep = '_')
      segms.df <- segms.df[segms.df$id %in% links.segments,]
    } else {
      warning("Object 'gaf.links' not defined, cannot filter based on link degree!!!")
    }
  }  
  
  ## Define segments ##
  #####################
  ## All segments
  nodes.df <- data.frame(x=c(rbind(segms.df$start, segms.df$end)),
                         y=0,
                         group=rep(1:nrow(segms.df), each=2),
                         rank=rep(segms.df$rank, each=2))
  segms.df$midpoint <- segms.df$start + ((segms.df$end - segms.df$start) / 2)
  ## Get haplotype segments
  if (!is.null(highlight.haplotype)) {
    hap.segms.df <- segms.df[segms.df$id %in% hap.segments,]
  }
  ## Define links ##
  ##################
  ## All links
  link.start <- ifelse(links$from.orient == '+', 
                       segms.df$end[match(links$from, segms.df$id)], 
                       segms.df$start[match(links$from, segms.df$id)])
  link.end <- ifelse(links$to.orient == '+', 
                     segms.df$start[match(links$to, segms.df$id)],
                     segms.df$end[match(links$to, segms.df$id)])
  #x.coords <- c(rbind(segms.df$end[match(links$from, segms.df$id)], segms.df$start[match(links$to, segms.df$id)]))
  x.coords <- c(rbind(link.start, link.end))
  y.coords <- c(rbind(segms.df$rank[match(links$from, segms.df$id)], segms.df$rank[match(links$to, segms.df$id)]))
  
  ## Impose arc height based on the segment order
  from <- match(links$from, segms.df$id)
  to <- match(links$to, segms.df$id)
  arc.height <- ( to - from ) - 1
 
  if (layout == 'linear') {
    arcs.df <- data.frame(x=rep(x.coords, each=2),
                          y=do.call(c, lapply(arc.height, function(x) c(0,x,x,0))),
                          group=rep(1:nrow(links), each=4))
  } else if (layout == 'offset') {
    arcs.df <- data.frame(x=rep(x.coords, each=2),
                          y=c(rbind(y.coords[c(TRUE, FALSE)], arc.height, arc.height, y.coords[c(FALSE, TRUE)])),
                          group=rep(1:nrow(links), each=4))
  }
  ## Add link IDs
  arcs.df$link.ids <- rep(links.ids, each=4)
  
  ## Add frequency to each link if defined and get a corrsponding link label
  if (!is.null(gaf.links)) {
    arcs.df$freq=rep(links$link.freq, each=4)
    ## Get link frequency labels
    link.label <- data.frame(xmin=segms.df$end[match(links$from, segms.df$id)], 
                              xmax=segms.df$start[match(links$to, segms.df$id)],
                              arc.height=arc.height,
                              freq=links$link.freq)
    link.label$midpoint <- link.label$xmin + ((link.label$xmax - link.label$xmin) / 2)
  }  
  
  ## Visualize nodes ##
  #####################
  if (layout == 'linear') {
    if (shape == 'rectangle') {
      segms.plt <- ggplot() +
        geom_rect(data=segms.df, aes(xmin=start, xmax=end, ymin=-0.4, ymax=0.4), size=2, colour = 'deepskyblue2', fill = 'deepskyblue2')
        #geom_text(data=segms.df, aes(x=midpoint, y=0.5, label=id), color='red')
    } else if (shape == 'roundrect') {
      segms.plt <- ggplot(nodes.df, aes(x = x, y = 0, group=group)) +
        geom_shape(radius = unit(0.25, 'cm'))
    }  
  } else if (layout == 'offset') {
    if (shape == 'rectangle') {
      segms.plt <- ggplot() +
        geom_rect(data=segms.df, aes(xmin=start, xmax=end, ymin=rank-0.4, ymax=rank + 0.4), size=2, colour = 'deepskyblue2', fill = 'deepskyblue2')
      #geom_text(data=segms.df, aes(x=midpoint, y=0.5, label=id), color='red')
    } else if (shape == 'roundrect') {
      segms.plt <- ggplot(nodes.df, aes(x = x, y = rank, group=group)) +
        geom_shape(radius = unit(0.25, 'cm'))
    }  
  }  
  
  ## Visualize links ##
  #####################
  if (is.null(link.frequency)) {
    final.plt <- segms.plt + 
      ggforce::geom_bezier(data=arcs.df, aes(x = x, y = y, group=group), arrow = arrow(type = arrow.head, length = unit(0.01, "npc")), inherit.aes = FALSE) +
      scale_x_continuous(labels = scales::comma)
  } else if (link.frequency == 'width') {  
    final.plt <- segms.plt +
      ggforce::geom_bezier(data=arcs.df, aes(x = x, y = y, group=group, size=freq), arrow = arrow(type = arrow.head, length = unit(0.01, "npc")), inherit.aes = FALSE) +
      geom_text(data=link.label, aes(x = midpoint, y = arc.height, label=freq), vjust=-0.2) +
      scale_x_continuous(labels = scales::comma) + scale_size_binned(range = c(0,2))
  } else if (link.frequency == 'color') {
    pal <- wesanderson::wes_palette(name = "Zissou1", n = max(link.label$freq), type = "continuous")
    final.plt <- segms.plt +
      ggforce::geom_bezier(data=arcs.df, aes(x = x, y = y, group=group, color=freq), arrow = arrow(type = arrow.head, length = unit(0.01, "npc")), inherit.aes = FALSE) +
      geom_text(data=link.label, aes(x = midpoint, y = arc.height, label=freq), vjust=-0.2) +
      scale_x_continuous(labels = scales::comma) + scale_size_binned(range = c(0,2)) +
      scale_color_gradientn(colours = pal)
  }  
				 
  ## Visualize Genes/Annotations ##
  #################################
  if (!is.null(gaf.annotation)) {
    ## Find where in gaps starts and stops of gene limits lie in our image
    gene_int <- findInterval(gaf.annotation$path.start,segms.df$start)
    gene_intEnd <- findInterval(gaf.annotation$path.end,segms.df$end)
    ## Space the starts and stops by # of gaps
    shifts <- gaf.annotation$path.start + ((gene_int - 1) * spacer)
    shiftsEnd <- gaf.annotation$path.end + ((gene_intEnd ) * spacer)
    gene_tbl <- unique(as.data.frame(cbind(as.numeric(shifts), as.numeric(shiftsEnd), gaf.annotation$q.name)))
    ## Create new dataframe for plotting
    gene_shift.gr <- GenomicRanges::GRanges(seqnames = 'genes', 
                                            ranges = IRanges::IRanges(start = as.numeric(gene_tbl$V1), end = as.numeric(gene_tbl$V2), id = gene_tbl$V3))
    gene_shift.gr$level <- disjointBins(gene_shift.gr) * -1
    gene.df <- as.data.frame(gene_shift.gr)
    gene.df$midpoint <- gene.df$start + ((gene.df$end - gene.df$start) / 2)
    final.plt <- final.plt +
      geom_rect(data=gene.df, aes(xmin=start, xmax=end, ymin=level-1, ymax=level-2), alpha=0.2, fill = 'darkblue') + 
      geom_text(data=gene.df, aes(x=midpoint, y=level-1.5, label = id), size=2, hjust = 0.5, color='black')
  }
  
  ## Visualize haplotype path ##
  ##############################
  if (!is.null(highlight.haplotype)) {
    if (is.null(link.frequency)) {
      final.plt <- final.plt +
        geom_rect(data=hap.segms.df, aes(xmin=start, xmax=end, ymin=rank-0.4, ymax=rank + 0.4), size=2, colour = 'red', fill = 'red') +
        ggforce::geom_bezier(data=arcs.df[arcs.df$link.ids %in% hap.links.ids,], aes(x = x, y = y, group=group), arrow = arrow(type = arrow.head, length = unit(0.01, "npc")), inherit.aes = FALSE, color='red')
    } else if (link.frequency == 'width') {  
      final.plt <- final.plt +
        geom_rect(data=hap.segms.df, aes(xmin=start, xmax=end, ymin=rank-0.4, ymax=rank + 0.4), size=2, colour = 'red', fill = 'red') +
        ggforce::geom_bezier(data=arcs.df[arcs.df$link.ids %in% hap.links.ids,], aes(x = x, y = y, group=group, size=freq), arrow = arrow(type = arrow.head, length = unit(0.01, "npc")), inherit.aes = FALSE, color='red')
    } else if (link.frequency == 'color') {
      pal <- wesanderson::wes_palette(name = "Zissou1", n = max(link.label$freq), type = "continuous")
      final.plt <- final.plt +
        geom_rect(data=hap.segms.df, aes(xmin=start, xmax=end, ymin=rank-0.4, ymax=rank + 0.4), size=2, colour = 'red', fill = 'red') +
        ggforce::geom_bezier(data=arcs.df[arcs.df$link.ids %in% hap.links.ids,], aes(x = x, y = y, group=group, color=freq), arrow = arrow(type = arrow.head, length = unit(0.01, "npc")), inherit.aes = FALSE, color='red')
    }  
  }  
  ## Apply theme
  graph.theme <- theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       axis.title.x=element_blank(),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       panel.background = element_blank())
  final.plt <- final.plt + graph.theme
  
  if (is.numeric(y.limit))
    final.plt <- final.plt + ylim(-y.limit,y.limit)
  
  ## Return final plotting object
  return(final.plt)
}
