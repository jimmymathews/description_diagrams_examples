# Contact: Jimmy  jimmy.mathews@alleninstitute.org
#
#' @import mclust

#' Title
#' 
#' Description
#' @param vertex_subset mm
#' @return blah blah
#' @export
calculate_minimal_diagram <- function(vertex_subset) {
}

calculate_all_minimal_diagrams <- function(dimension = 3) {
	number_subsets <- 2^(2^dimension)
	subsets <- t(sapply(c(1:number_subsets), FUN=function(s) { return(as.integer(intToBits(s))[1:(2^dimension)]) }))
	return(subsets)
}

