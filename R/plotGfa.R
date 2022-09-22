#' Plot GFA from loaded data table
#' 
#' This function takes an output from readGfa.R function and makes a graph plot (Adjust annotation!!!)
#'
#' @param gfa.tbl A path to a GFA file containing sequence graph.
#' @param min.segment.length A minimum size of segment to be plotted.
#' @param spacer.width User defined fraction to the total segment length to be used as node spacer.
#' @param order.by Define a column to be used for node ordering. [TODO]
#' @param layout Overall layout of the graph, either 'linear' or 'circular'.
#' @param shape Either 'rectangle' or 'roundrect' to plot graph nodes.
#' @param arrow.head Default 'closed' TODO what are the other options?
#' @param gaf.links A \code{tibble} table containing links present in each haplotype after alignment to the graph.
#' @param link.frequency Visualize link frequency either by the thickness of a link ('width') or a color gradient ('color').
#' @importFrom GenomicRanges shift GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors lapply
#' @importFrom ggforce geom_bezier
#' @author David Porubsky, Sean McGee & Karynne Patterson
#' @export
#' 
plotGfa <- function(gfa.tbl=NULL, min.segment.length=0, spacer.width=0.05, order.by='offset', layout='linear', shape='rectangle', arrow.head='closed', gaf.links=NULL, link.frequency=NULL) {
  ## Check user input ##
  ######################
  ## Get data from loaded GFA file
  segments <- gfa.tbl$segments
  links <- gfa.tbl$links
  
  ## Get link frequency from links if gaf.links are defined
  if (!is.null(gaf.links)) {
    gaf.links.ids <- paste(gaf.links$from, gaf.links$to, sep = '_')
    link.freq <- as.data.frame(table(gaf.links.ids))
    link.freq <- link.freq[order(link.freq$Freq, decreasing = TRUE),]
  }
  links.ids <- paste(links$from, links$to, sep = '_')
  links$link.freq <- link.freq$Freq[match(links.ids, link.freq$gaf.links.ids)]
  
  ## Filter data ##
  #################
  ## Filter segments by size [TODO]
  #segments <- segments[segments$LN >= min.segment.length,]
  #links <- links[links$from %in% segments$segment.id & links$to %in% segments$segment.id,]
  
  ### Prepare data for plotting ###
  #################################
  ## Define segment spacer as fraction of the total lenght of all segments
  spacer <- sum(segments$LN) * spacer.width
  
  ## Order segments ##
  if (order.by == 'offset') {
    segments <- segments[order(segments$SO),]
  }  
  
  ## Space out graph segments into a single line ##
  segms.gr <- GenomicRanges::GRanges(seqnames = 'nodes', 
                                     ranges = IRanges::IRanges(start = 1, end = segments$LN), id = segments$segment.id, rank=segments$SR)
  shifts <- GenomicRanges::width(segms.gr)
  shifts <- cumsum(shifts + spacer)
  segms.gr[-1] <- GenomicRanges::shift(segms.gr[-1], shift = shifts[-length(shifts)])
  
  ## Define segments ##
  #####################
  segms.df <- as.data.frame(segms.gr)
  nodes.df <- data.frame(x=c(rbind(segms.df$start, segms.df$end)),
                         y=0,
                         group=rep(1:nrow(segms.df), each=2),
                         rank=rep(segms.df$rank, each=2))
  segms.df$midpoint <- segms.df$start + ((segms.df$end - segms.df$start) / 2)
  
  ## Define links ##
  ##################
  link.start <- ifelse(links$from.orient == '+', 
                       segms.df$end[match(links$from, segms.df$id)], 
                       segms.df$start[match(links$from, segms.df$id)])
  link.end <- ifelse(links$to.orient == '+', 
                     segms.df$start[match(links$to, segms.df$id)],
                     segms.df$end[match(links$to, segms.df$id)])
  x.coords <- c(rbind(segms.df$end[match(links$from, segms.df$id)], segms.df$start[match(links$to, segms.df$id)]))
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
  
  ## Apply theme
  graph.theme <- theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       axis.title.x=element_blank(),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       panel.background = element_blank())
  final.plt <- final.plt + graph.theme
  
  ## Return final plotting object
  return(final.plt)
}