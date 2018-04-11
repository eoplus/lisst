
#' \pkg{lisst}: Utilities for analysis of LISST data
#'
#' An R package to read, manipulate and visualize data from the Laser In-Situ 
#' Scattering and Transmissometry (LISST) instruments, Sequoia Scientific, Inc.
#'
#' @author Alexandre Castagna
#'
#' @docType package
#' @import units
#' @import errors
#' @import quantities
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
	inherits(x, "lisst") && any(!sapply(.lattributes(x), is.null) && !is.null(x$Time))
}

#' drop lisst
#'
#' Drop lisst class and attributes.
#'
#' @param x A lisst object.
#'
#' @return A \code{data.frame}, with all columns of class \code{quantities} in 
#' exception of those that contain time information.
#'
#' @examples
#' x <- drop_lisst(donkmeer_bin)
#' attributes(x)
#'
#' @export

drop_lisst <- function(x) {
	structure(x, type = NULL, lproc = NULL, linst = NULL, lmodl = NULL, 
		zscat = NULL, lxts = NULL, class = "data.frame")
}


