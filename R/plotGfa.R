
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
#' @importFrom GenomicRanges shift GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors lapply
#' @importFrom ggforce geom_bezier
#' @author David Porubsky
#' @export
#' 
plotGfa <- function(gfa.tbl=NULL, min.segment.length=0, spacer.width=0.05, order.by='offset', layout='linear', shape='rectangle') {
  ## Check user input ##
  segments <- gfa.tbl$segments
  links <- gfa.tbl$links
  
  ## Filter data ##
  ## Filter segments by size
  segments <- segments[segments$LN >= min.segment.length,]
  links <- links[links$from %in% segments$segment.id & links$to %in% segments$segment.id,]
  
  ## Define segment spacer as fraction of the total lenght of all segments
  spacer.width <- sum(segments$LN) * 0.05
  
  ## Order segments ##
  if (order.by == 'offset') {
    segments <- segments[order(segments$SO),]
  }  
  
  ## Space out graph segments into a single line ##
  segms.gr <- GenomicRanges::GRanges(seqnames = 'nodes', 
                                     ranges = IRanges::IRanges(start = 1, end = segments$LN), id = segments$segment.id, rank=segments$SR)
  shifts <- GenomicRanges::width(segms.gr)
  shifts <- cumsum(shifts + spacer.width)
  segms.gr[-1] <- GenomicRanges::shift(segms.gr[-1], shift = shifts[-length(shifts)])
  
  ## Prepare data for plotting ##
  ## Define segments
  segms.df <- as.data.frame(segms.gr)
  nodes.df <- data.frame(x=c(rbind(segms.df$start, segms.df$end)),
                         y=0,
                         group=rep(1:nrow(segms.df), each=2),
                         rank=rep(segms.df$rank, each=2))
  segms.df$midpoint <- segms.df$start + ((segms.df$end - segms.df$start) / 2)
  
  ## Define links ##
  x.coords <- c(rbind(segms.df$end[match(links$from, segms.df$id)], segms.df$start[match(links$to, segms.df$id)]))
  y.coords <- c(rbind(segms.df$rank[match(links$from, segms.df$id)], segms.df$rank[match(links$to, segms.df$id)]))
  
  from <- as.numeric(gsub(links$from, pattern = 's', replacement = ''))
  to <- as.numeric(gsub(links$to, pattern = 's', replacement = ''))
  arc.height <- ( to - from ) - 1
  if (layout == 'linear') {
    arcs.df <- data.frame(x=rep(x.coords, each=2),
                          y=do.call(c, lapply(arc.height, function(x) c(0,x,x,0))),
                          group=rep(1:nrow(links), each=4))
  } else if (layout == 'offset') {
    arcs.df <- data.frame(x=rep(x.coords, each=2),
                          y=c(rbind(segms.df$rank[links$from], arc.height, arc.height, segms.df$rank[links$to])),
                          group=rep(1:nrow(links), each=4))
  }  
  
  ## Visualize nodes ##
  if (layout == 'linear') {
    if (shape == 'rectangle') {
      segms.plt <- ggplot() +
        geom_rect(data=segms.df, aes(xmin=start, xmax=end, ymin=-0.4, ymax=0.4), colour = 'black', fill = 'black')
        #geom_text(data=segms.df, aes(x=midpoint, y=0.5, label=id), color='red')
    } else if (shape == 'roundrect') {
      segms.plt <- ggplot(nodes.df, aes(x = x, y = 0, group=group)) +
        geom_shape(radius = unit(0.25, 'cm'))
    }  
  } else if (layout == 'offset') {
    if (shape == 'rectangle') {
      segms.plt <- ggplot() +
        geom_rect(data=segms.df, aes(xmin=start, xmax=end, ymin=rank-0.4, ymax=rank + 0.4), colour = 'black', fill = 'black')
      #geom_text(data=segms.df, aes(x=midpoint, y=0.5, label=id), color='red')
    } else if (shape == 'roundrect') {
      segms.plt <- ggplot(nodes.df, aes(x = x, y = rank, group=group)) +
        geom_shape(radius = unit(0.25, 'cm'))
    }  
  }  
  
  ## Visualize links ##
  final.plt <- segms.plt + 
    ggforce::geom_bezier(data=arcs.df, aes(x = x, y = y, group=group), arrow = arrow(length = unit(0.01, "npc")), inherit.aes = FALSE)
  
  ## Return final plotting object
  return(final.plt)
}
