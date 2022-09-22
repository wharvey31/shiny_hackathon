# returns the max absolute value of a collection of primitive numerics 
max_absolute_value <- function(primitive_numerics) {
	max_abs_value = abs(primitive_numerics[1])
	for (x in 2:length(primitive_numerics)) {
		abs_primitive_numeric <- abs(primitive_numerics[x])
		if (abs_primitive_numeric > max_abs_value)
			max_abs_value <- abs_primitive_numeric
	}
	return(max_abs_value)
}

# 
segments_for_haplotype_links <- function(all_segments, haplotype_links){
	haplotype_segment_ids <- (c(haplotype_links$from, haplotype_links$to))
	haplotype_segment_ids_set <- haplotype_segment_ids[!duplicated(haplotype_segment_ids)]
	segments_for_haplotype_links <- subset(all_segments, (segment.id %in% haplotype_segment_ids_set))
	return(segments_for_haplotype_links)
}

rgb_to_hex <- function(rgb_comm){
  rgb_comm = strsplit(split = ",", x = rgb_comm) %>% unlist()
  return(rgb(red = rgb_comm[1], green = rgb_comm[2], blue = rgb_comm[3], maxColorValue = 255))
}

