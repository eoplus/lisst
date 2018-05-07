
#' Fit a PSD model to the lisst data
#'
#' The function performs a fit of the available models of PSD.
#' 
#' @param x     A PSD ('vol' or 'pnc') lisst object.
#' @param model Curently ignored. Only Junge model implemented.
#'
#' @export

lfit <- function(x, model) {
	stopifnot(is.lisst(x))
	typ <- attr(x, 'type')
	if(!(typ == 'vol' || typ == 'pnc'))
		stop("A PSD model fit can only be performed in a 'vol' or 'pnc' lisst object", 
			call. = FALSE)
	x <- lget(x, 'pnc')

	linst <- attr(x, "linst")
	lmodl <- attr(x, "lmodl")
	lproc <- attr(x, "lproc")

	bins  <- as.numeric(lmodl$binr[[lproc$ity]][, 3])
	junge <- numeric(nrow(x))
	mat   <- as.matrix(x[, 1:lmodl$nring])
	id    <- which(mat[, 1] != 0)
	for(i in id) {
		jfit <- lm(log10(mat[i, ])~log10(bins))
		junge[i] <- as.vector(coefficients(jfit)[2])
	}
	return(junge)
}

#' Descriptive statistics for lisst objects
#'
#' The functions calculate the average or the median per variable of a lisst 
#' object, for the whole data or in intervals. The dispersion of the dataset or
#' intervals is also retrieved and stored.
#'
#' @param x    A lisst object.
#' @param brks A vector with the breaks (intervals) for the aggregation. See 
#'             details.
#' @param fun  A function to perform the aggregation (mean, median, sd). 
#'             Defaults to mean.
#' @param ...  Arguments to be passed to aggregation functions.
#'
#' @details 
#' The breaks (intervals) are passed directly to the subset function, so must 
#' now be supplied in final form.
#'
#' @export

lstat <- function(x, brks, fun = 'mean', ...) {
	stopifnot(is.lisst(x))
	if(missing(brks)) stop("brks missing, with no default", call. = FALSE)
	if(missing(fun)) fun <- 'mean'
	else if(!(fun == 'mean' || fun == 'median' || fun == 'sd'))
		stop("fun must be 'mean', 'median' or 'sd'.", call. = FALSE)

	xl <- list()
	for(i in 1:length(brks)) xl[[i]] <- x[brks[[i]], ]
	xb <- do.call(rbind, lapply(xl, fun, ...))
	rownames(xb) <- 1:nrow(xb)
	xb
}

#' @describeIn lstat Compute the mean for lisst objects
#'
#' @examples
#' mean(donkmeer_bin)
#'
#' @export

mean.lisst <- function(x, ...) {
	stopifnot(is.lisst(x))
	xm <- x[1, , drop = FALSE]
	xm[1, ] <- as.data.frame(sapply(x, mean, na.rm = TRUE, simplify = F))
	xm
}

#' @describeIn lstat Compute the median for lisst objects
#'
#' @examples
#' median(donkmeer_bin)
#'
#' @export

median.lisst <- function(x, ...) {
	stopifnot(is.lisst(x))
	xm <- x[1, , drop = FALSE]
	xm[1, ] <- as.data.frame(sapply(x, median, na.rm = TRUE, simplify = F))
	xm
}

#' An S3 generic for sd
#'@export

sd <- function(x, ...) UseMethod("sd")
sd.default <- stats::sd

#' @describeIn lstat Compute the median for lisst objects
#'
#' @examples
#' median(donkmeer_bin)
#'
#' @export

sd.lisst <- function(x, ...) {
	stopifnot(is.lisst(x))
	xm <- x[1, , drop = FALSE]
	xm$Time <- numeric(1)
	xm[1, ] <- as.data.frame(sapply(x, sd, na.rm = TRUE, simplify = F))
	xm
}

