library(ggplot2)
library(gggenes)
library(tidyverse)
library(glue)


#'is_hex = function to check if a column is hex
is_hex = function(color_str){
  return( grepl(pattern =  "^#[0-9A-F]{6}$", x = color_str ) )
}

#'is_rgb = function to check if a column is in rgb
is_rgb = function(color_str){
  return(strsplit(color_str, split = ",") %>% unlist() %>% grepl(pattern = "^[0-9]") %>% sum() == 3 )
}

#'rgb_to_hex : convert comma_sep rgb string to hex color.
rgb_to_hex = function(rgb_comm){
  rgb_comm = strsplit(split = ",", x = rgb_comm) %>% unlist()
  return(rgb(red = rgb_comm[1], green = rgb_comm[2], blue = rgb_comm[3], maxColorValue = 255))
}

#' check_bed : check that certain columns are as expected for a bed.
check_bed = function(bed_df){
  #check if empty
  if(any( dim(bed_df) == 0 )){
    stop(glue("Bed file is empty.") )
  }
  #check if not enough columns
  stopifnot( "Not enough columns passed in bed. Require 3 or more" = ncol(bed_df) >= 3)
  #check that start and stop are integers
  stopifnot( "start and stop columns are not integers." = 
               all( bed_df[,c(2,3), drop = T] %>% summarize_all(class) %>% unlist() == "numeric" ))
  #check that start < stop
  stopifnot("start column is not always less than or equal to stop column" = 
              sum(bed_df[,2] <= bed_df[,3]) == nrow(bed_df)
              )
}

#'load_bed : load a bed for annotations and make columns needed for plot_bed_annot_track func.
load_annotation_bed = function(bed_path, color_col = NULL){
  STRAND_COL = 6 #column of strand in bed file
  out_df = read_tsv(file = bed_path , col_names = F, comment = "#")
  check_bed(out_df)
  colnames(out_df)[1:3] = c("contig", "start", "stop") 
  #add direction:
  if(ncol(out_df) >= STRAND_COL){
    if( all( levels(as.factor(out_df[,STRAND_COL, drop = T])) %in% c( "+", "-", "<", ">", "forward", "reverse", "for", "rev" ) ) ){
      colnames(out_df)[STRAND_COL] = 'strand' 
      out_df = out_df %>% mutate( strand = ifelse(strand %in% c("+", ">", "forward", "for"), yes = "+", no = "-"))
    }
  }
  #add color
  if( !(is.null(color_col)) ){
    if(all(sapply(out_df[,color_col, drop = T], FUN = is_rgb))){ #check if they're RGB
      out_df$color = sapply(out_df[,color_col, drop = T], rgb_to_hex)
      names(out_df$color) = NULL
    }else{
      out_df$color = out_df[,color_col, drop = T]
    }
    
  }
  stopifnot(is.character(out_df$contig), is.numeric(out_df$start), is.numeric(out_df$stop))
  return(out_df)
}


#'plot_bed_annot_track : generate a gggenes plot of a genomic annotation, or add a track to an existing plot
#'input_mandatory: bed dataframe with minimum columns: contig , start, stop
#'input_mandatory: track_name: name you want to give track
#'input_optional: p : existing plot to add track to.
#'output: ggplot p
#'NOTE: if want to add colors, include column called color in bed_df
#'NOTE: assuming facet wrapping is by contig.
plot_bed_annot_track = function(track_name, bed_df, p = NULL, facet_col = "contig"){
  if(is.null(p)){
    p = ggplot()
  }
  if("color" %in% colnames(bed_df)){
    if(all(sapply(bed_df$color, FUN = is_hex) )){
      fill_col = bed_df$color
    }else{
      fill_col = "black"
      #fill_col = as.factor(bed_df$color)
    }
  }else{
    fill_col = "red"
  }
  
  p = p + gggenes::geom_gene_arrow(data = bed_df , mapping =  aes(xmin = start, xmax = stop, y = track_name), 
                                   fill = fill_col) +
    theme(legend.position="none")
  
  #check if multiple contigs
  if( length(unique(bed_df$contig)) > 1 ){
    p = p + facet_wrap(~contig, ncol = 1) + theme(legend.position="none")
  }
  #
  p
}

#test
# dup_bed1 = load_annotation_bed(bed_path = "example_data/HG01071_2/HG01071_2.duplicons.bed", color_col = 9)
# trf_bed1 = load_annotation_bed(bed_path = "example_data/HG01071_2/HG01071_2.trf.bed")
# rm_bed1 = load_annotation_bed(bed_path = "example_data/HG01071_2/HG01071_2.rm.bed", color_col = 4)
# 
# 
# p = plot_bed_annot_track(track_name = "dups", bed_df = dup_bed1)
# p = plot_bed_annot_track(track_name = "trf", bed_df = trf_bed1 , p = p)
# p = plot_bed_annot_track(track_name = "repeat_masker", bed_df = rm_bed1, p  = p)
# 
# 
# # #try facet
# dup_bed1 = load_annotation_bed(bed_path = "example_data/combined/combined_duplicons.bed", color_col = 9)
# trf_bed1 = load_annotation_bed(bed_path = "example_data/combined/combined_trf.bed")
# rm_bed1 = load_annotation_bed(bed_path = "example_data/combined/combined_rm.bed", color_col = 4)
# 
# p = plot_bed_annot_track(track_name = "dups", bed_df = dup_bed1)
# p = plot_bed_annot_track(track_name = "trf", bed_df = trf_bed1 , p = p)
# p = plot_bed_annot_track(track_name = "repeat_masker", bed_df = rm_bed1, p  = p)
