#' Add plotting coordinates to each segment taken from GFA file
#' 
#' This function takes an output from readGfa.R function and makes and adds segment plotting coordinates to the final data object.
#'
#' @inheritParams plotGfa
#' @importFrom GenomicRanges shift GRanges
#' @importFrom IRanges IRanges
#' @author David Porubsky
#' @export
#' 
addGraphPlottingCoords <- function(gfa.tbl=NULL, spacer.width=0.05) {
  ## Check user input ##
  ######################
  ## Get data from loaded GFA file
  segments <- gfa.tbl$segments
  
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
  segms.df <- as.data.frame(segms.gr)
  ## Add plotting coordinates to GFA data object
  gfa.tbl$graph.plt.coords <-  segms.df
  return(gfa.tbl)
}
