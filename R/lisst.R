
#' \pkg{lisst}: Utilities for analysis of LISST data
#'
#' An R package to read, manipulate and visualize data from the Laser In-Situ 
#' Scattering and Transmissometry (LISST) instruments, Sequoia Scientific, Inc.
#'
#' @author Alexandre Castagna
#'
#' @docType package
#' @import stats
#' @import units
#' @import xts
#' @name lisst-package
NULL

#' Test for a lisst object
#'
#' Check that an object loosely conforms to a lisst object.
#'
#' @param x An object to be tested.
#'
#' @return A logical value indicating if the argument is a lisst object.
#'
#' @examples
#' is.lisst(donkmeer_pro)
#' is(donkmeer_pro, 'lisst')
#'
#' @export

is.lisst <- function(x) {
	inherits(x, "lisst") && any(!sapply(.lattributes(x), is.null))
}

#' drop lisst
#'
#' Drop lisst class and attributes.
#'
#' @param x A lisst object.
#'
#' @return A \code{data.frame}, with all columns of class \code{units} in 
#' exception of those that contain time information.
#'
#' @examples
#' x <- drop_lisst(donkmeer_bin)
#' attributes(x)
#'
#' @export

drop_lisst <- function(x) {
	rownames(x) <- 1:nrow(x)
	structure(x, type = NULL, lproc = NULL, linst = NULL, lmodl = NULL, 
		zscat = NULL, lxts = NULL, class = "data.frame")
}

#' Extract size range of bins
#'
#' Extracts the size range and median size per bin per LISST model and particle
#' inversion type.
#'
#' @param x   A lisst object.
#' @param mod A LISST model ('100(X)B', '100(X)C' or '200X').
#' @param ity Inversion type ('ss' or 'rs').
#'
#' @details
#' If x is provided, mod and ity (if present) will be extracted from the 
#' object's metadata. If provided, mod and ity will take precedence. Note that
#' since inversion is not implemented in this version, only lisst objects 
#' created from LISST SOP processed files will contain ity information.
#'
#' @return A \code{units} \code{matrix} with lower and upper size limits and 
#' median size in µm.
#'
#' @examples
#' lbinr(mooring)
#' lbinr(,'100XC', 'rs')
#'
#' @export

lbinr <- function(x, ity) {
	if(is(x, 'lisst')) {
		lmodl <- attr(x, 'lmodl')
	} else if(is(x, 'character')) {
		x     <- sub("X", "", x, ignore.case = TRUE)
		lmodl <- switch(x,
			"200"  = .getmodp(list(mod = "200", dty = "A")),
			"100B" = .getmodp(list(mod = "100", dty = "B")),
			"100C" = .getmodp(list(mod = "100", dty = "C")),
			stop("model must be one of '100(X)B', '100(X)C' or ", 
				"'200X'", call. = FALSE)
		)
	} else {
		stop('x must be a lisst or a character object', call. = FALSE)
	}
	if(missing(ity)) {
		ity <- attr(x, 'lproc')$ity
		if(ifelse(is.null(ity), TRUE, is.na(ity))) {
			if(lmodl$mod == '200')
				ity <- 'ss'
			else
				stop('ity is not defined in the lisst object ',
					'or needs more specification than just model', 
					call. = FALSE)
		}
	}
	if(is.na(ity) || !(ity == 'ss' || ity == 'rs')) {
		stop("ity must be either 'ss' or 'rs'", call. = FALSE)
	} else {
		if(lmodl$mod == '100') cat(paste('ity:', ity), '\n') 
		lmodl$binr[[ity]]
	}
}

#' Extract angle range of bins
#'
#' Extracts the (in water) angle range and median angle per bin per LISST model.
#'
#' @param x A lisst object or a LISST model ('100(X)B', '100(X)C' or '200X').
#'
#' @return A \code{units} \code{matrix} with lower and upper angle limits and 
#' median angle in radians.
#'
#' @examples
#' lwang(mooring)
#' lwang('100XC')
#'
#' @export

lwang <- function(x) {
	if(is(x, 'lisst')) {
		lmodl <- attr(x, 'lmodl')
	} else if(is(x, 'character')) {
		x     <- sub("X", "", x, ignore.case = TRUE)
		lmodl <- switch(x,
			"200"  = .getmodp(list(mod = "200", dty = "A")),
			"100B" = .getmodp(list(mod = "100", dty = "B")),
			"100C" = .getmodp(list(mod = "100", dty = "C")),
			stop("model must be one of '100(X)B', '100(X)C' or ", 
				"'200X'", call. = FALSE)
		)
	} else {
		stop('x must be a lisst or a character object')
	}
	lmodl$wang
}

