#' Read GAF from an input file
#' 
#' This function takes an GFA output file from minigraph and loads the file
#' into a set of segments and links.
#'
#' @param gfa.file A path to a GAF file containing sequence graph.
#' @param store.fasta Set to \code{TRUE} if the FASTA sequence should be retained and exported.
#' @param restrict.gfa.tags Define a set of GAF tag ids (e.g. LN) to be reported in the output. 
#' @importFrom stringr str_split
#' @importFrom dplyr bind_cols bind_rows
#' @importFrom S4Vectors lapply
#' @return A \code{list} object containing separate table for segments, links and/or FASTA sequences if 'store.fasta' is set to \code{TRUE}
#' @author David Porubsky
#' @export
#' 
readGfa <- function(gfa.file, store.fasta=FALSE, restrict.gfa.tags=NULL) {

  ## Check user input ##
  if (is.null(gfa.file)) {
    stop("Path to a PAF file to load is not defined!!!")
  }

  ## Check if file exists and is TAB delimited 
  if (file.exists(gfa.file)) {
    con <- file(gfa.file,"r")
    first.line <- readLines(con, n=1)
    n.fields <- length(stringr::str_split(first.line, "\t")[[1]])
    if (n.fields < 3) {
      stop("User defined 'gfa.file' has less then 3 expected tab-delimeted fields!!!")
    }
    close(con)
  } else {
    stop("User defined 'gfa.file' doesn't exist!!!")
  }
  
  ## Load GAF file ##
  file.con <- file(gfa.file, "r")
  
  segments <- list()
  links <- list()
  fastas <- list()
  while (TRUE) {
    line <- readLines(file.con, n = 1)
    if (length(line) != 0) {
      fields <- stringr::str_split(line, pattern = "\t")
      ## Get gfa record.type
      record.type <-  fields[[1]][1]
      if (record.type == 'S') {
        segm.fields <- S4Vectors::lapply(fields, "[", c(1:2))
        field.names <-  c('record.type', 'segment.id')
        attr(segm.fields[[1]], "names") <- field.names
        ## Parse GFA tags
        gfa.tags <- S4Vectors::lapply(fields, function(x) paste(x[4:length(x)]))
        ## Merge objects
        segms <- dplyr::bind_cols(dplyr::bind_rows(segm.fields[[1]]), processGfaTags(gfa.tags = gfa.tags))
        segments[[length(segments) + 1]] <- segms
        ## Keep FASTA if desired
        if (store.fasta) {
          fa.obj <- Biostrings::DNAStringSet(fields[[1]][3])
          names(fa.obj) <- segms$segment.id
          fastas[[length(fastas) + 1]] <- fa.obj
        }
      } else if (record.type == 'L') {
        link.fields <- S4Vectors::lapply(fields, "[", c(1:6))
        field.names <-  c('record.type', 'from', 'from.orient', 'to', 'to.orient', 'cg.overlap')
        attr(link.fields[[1]], "names") <- field.names
        ## Parse GFA tags
        gfa.tags <- S4Vectors::lapply(fields, function(x) paste(x[7:length(x)]))
        ## Merge objects
        lnks <- dplyr::bind_cols(dplyr::bind_rows(link.fields[[1]]), processGfaTags(gfa.tags = gfa.tags, restrict.gfa.tags = restrict.gfa.tags))
        links[[length(links) + 1]] <- lnks
      }
    } else {
      break
    }
  }  
  close(file.con)
  segments.tbl <- dplyr::bind_rows(segments)
  links.tbl <- dplyr::bind_rows(links)
  fastas <- do.call(c, fastas)
  ## Return list of loaded data objects
  if (store.fasta) {
    list(segments = segments.tbl, links = links.tbl, fastas = fastas)
  } else {
    list(segments = segments.tbl, links = links.tbl)
  }  
}


#' Process GFA specific alignment tags.
#'
#' @param pfa.tags A \code{list} of GFA specific tags.
#' @inheritParams readGaf
#' @importFrom stringr str_split
#' @importFrom dplyr bind_rows
#' @author David Porubsky
#' @export
processGfaTags <- function(gfa.tags, restrict.gfa.tags=NULL) {
  ## Takes only expected tag values
  allowed.tags <- c('LN', 'SN', 'SO', 'SR', 'L1', 'L2')
  if (!is.null(restrict.gfa.tags)) {
    restrict.gfa.tags <- restrict.gfa.tags[restrict.gfa.tags %in% allowed.tags]
    if (length(restrict.gfa.tags) == 0) {
      message(paste0("Submitted 'restrict.paf.tags' are not present in the allowed set of PAF tags: ", paste(allowed.tags, collapse = '; ')))
      restrict.gfa.tags <- allowed.tags
    }
  } else {
    restrict.gfa.tags <- allowed.tags 
  }  
  
  tags <- character()
  to.numeric <- integer()
  res <- list()
  n <- length(gfa.tags)
  t.idx <- 0

  split.tags <- stringr::str_split(gfa.tags[[1]], ":")
  for (tag in split.tags) {
    if (!is.null(restrict.gfa.tags) & tag[1] %in% restrict.gfa.tags) {
      if ( !(tag[1] %in% tags)) {
        t.idx <- t.idx + 1
        if ( tag[2] %in% c("f", "H", "i")) {
          to.numeric <- c(to.numeric, t.idx)
        }
        res[[ tag[1] ]] <- rep(NA, n)
        tags <- c(tags, tag[1])
      }
      res[[ tag[1] ]][1] <- tag[3]
    }  
  }
  
  for(i in to.numeric){
    res[[i]] <- as.numeric(res[[i]])
  }
  dplyr::bind_rows(res)
}
