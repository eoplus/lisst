# TODO
# add an subset and apply function that preserves structure, units and attributes.

#' Fit a PSD model to the lisst data
#'
#' The function performs a fit of the available models of PSD.
#' 
#' @param lo   A PSD ('vol' or 'pnc') lisst object.
#' @param type Curently ignored.
#'
#' @export


lgetfit <- function(lo, type) {
	if(!is(lo, "lisst"))
		stop("lo must be a lisst object", call. = FALSE)
	typ <- attr(lo, "type")
	if(!(typ == "vol" || typ == "pnc"))
		stop("lo must be a PSD ('vol' or 'pnc') lisst object", call. = FALSE)

	if(typ == "vol") lo <- lgetpnc(lo)

	linst <- attr(lo, "linst")
	lmodl <- attr(lo, "lmodl")

	bins  <- as.numeric(lmodl$binr[[linst$ity]][, 3])
	junge <- numeric(nrow(lo))
	mat   <- as.matrix(lo[, 1:lmodl$nring])
	id    <- which(mat[, 1] != 0)
	for(i in id) {
		jfit <- lm(log10(mat[i, ])~log10(bins))
		junge[i] <- as.vector(coefficients(jfit)[2])
	}
	return(junge)
}


#
# Function: getnconc
#
