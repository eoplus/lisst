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


# If any kind of mean is made, is important to provide back the n in each class and the variability.
# should keep breaks as well...

lgetstat <- function(lo, Depth, Time) {
	if(!missing(Depth) && !missing(Time))
		stop("Depth and Time cannot be specified at the same time", call. = FALSE)
	if(missing(Depth) && missing(Time)) Depth <- range(lo$Depth)
	else if(!is.numeric(Depth))
		stop("depth must be numeric", call. = FALSE)
	else if(length(Depth) < 2)
		stop("depth must have length greater than 1", call. = FALSE)
	else if(any(is.na(Depth)))
		stop("NAs not allowed in depth argument", call. = FALSE)
	else
		units(Depth) <- "m"
	param <- "Depth"

	if(!missing(Time)) {
		if(length(Time) < 2)
			stop("time must have length greater than 1", call. = FALSE)
		else if(any(is.na(Time)))
			stop("NAs not allowed in time argument", call. = FALSE)
		else if(!is(Time, "POSIXct")) Time <- as.POSIXct(Time)
		param <- "Time"
	}

	idx <- numeric(nrow(lo))
	for(i in 1:length(Depth)) {
		idx <- idx + (lo[, param] > get(param)[i])
	}
	lom <- los <- lo[1:max(idx), ]
	for(i in 1:max(idx)) {
		id <- idx == i
		lom[i, ] <- lgetmean(subset(lo, id))
		los[i, ] <- lgetsd(subset(lo, id))
	}
	blon <- by(lo[id0, ], idx[id0], nrow)

	lmodl <- attr(lo, "lmodl")
	lom <- lo[1:length(blom), , drop = FALSE]
	los <- matrix(NA, ncol = ncol(lo), nrow = length(blom))
	colnames(los) <- c(lmodl$lvarn, "Time")
	for(i in 1:length(bylom)) {	
		lom[i, ] <- bylom[[i]]
		los[i, ] <- t(bylos[[i]])
	}
	attr(lom, "lproc") <- c(attr(lom, "lproc"), list(sd = los, n = as.vector(bylon)))
	lom
}


lgetmean <- function(lo) {
	lom <- lo[1, , drop = FALSE]
	lom[1, ] <- as.data.frame(sapply(lo, mean, simplify = F))
	lom
}

lgetsd <- function(lo) {
	unlist(sapply(lo, sd, simplify = F))
}

#
# Function: getnconc
#
