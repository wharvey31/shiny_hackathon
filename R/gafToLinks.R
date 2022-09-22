#' Load GFA file into a set of graph links
#'
#' @param gfa.file A path to a GFA file containing sequence graph.
#' @importFrom GenomicRanges GRanges sort
#' @importFrom IRanges IRanges
#' @importFrom dplyr bind_cols bind_rows
#' @author David Porubsky
#' @export
#' 
gafToLinks <- function(gaf.file=NULL) {
  ## Check user input ##
  if (is.null(gaf.file)) {
    stop("Path to a GAF file to load is not defined!!!")
  }
  
  ## Read GAF file
  gaf.df <- readGaf(gaf.file = gaf.file)
  ## Remove NAs
  #gaf.df <- gaf.df[!is.na(gaf.df$s.start),]
  
  ## Covert to Genomic ranges
  #gaf.gr <- GenomicRanges::GRanges(seqnames = gaf.df$q.name, ranges=IRanges::IRanges(start=gaf.df$p.start, end=gaf.df$p.end), strand=gaf.df$s.strand, id=gaf.df$s.name)
  #gaf.gr <- GenomicRanges::sort(gaf.gr)
  #gaf.grl <- split(gaf.gr, seqnames(gaf.gr))
  gaf.l <- split(gaf.df, paste(gaf.df$q.name, gaf.df$record, sep = '_'))
  
  ## Get subsequent links
  all.links <- list()
  for (i in seq_along(gaf.l)) {
    path <- gaf.l[[i]]
    seq.name <- names(gaf.l[i])
    seg.ids <- path$s.name
    ## Get all subsequent pairs
    from <- seg.ids[-length(seg.ids)]
    to <- seg.ids[-1]
    from.orient <-path$s.strand[-length(seg.ids)]
    to.orient <- path$s.strand[-1]
    link.ids <- paste(from, to, sep = '_')
    ## Construct link table
    links <- dplyr::bind_cols(record.type='L', from=from, from.orient=from.orient, to=to, to.orient=to.orient, SN=seq.name)
    all.links[[i]] <-  links
  }
  all.links <-  dplyr::bind_rows(all.links)

  # Add halpotype
  haplotype_name <- function(row, output) {
	  underscore_locations <- unlist(gregexpr("_", row[6]))
	  haplotype_name <- substr(row[6], 1, underscore_locations[2] -1)
	  return(haplotype_name)
  }
  all.links$haplotype <- apply(all.links, 1, haplotype_name) 
  
  
  ## Return links object
  return(all.links)
}
