

#
# Function: "["
#

`[.lisst` <- function(x, i, j, ..., drop = TRUE) {
	if(((missing(i) && length(j) == 1) || (missing(j) && length(i) == 1)) && drop) {
		NextMethod(drop = drop)
	} else {
		structure(NextMethod(drop = drop), 
			"sn" = attr(x, "sn"), 
			"type" = attr(x, "type"),
			"linst" = attr(x, "linst"),
			"lmodl" = attr(x, "lmodl"),
			"zscat" = attr(x, "zscat"),
			class = c("lisst", "data.frame")
		)
	}
}

