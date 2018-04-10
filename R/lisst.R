
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
#' @name lisst-package
NULL

#' drop lisst
#'
#' Drop lisst class and attributes.
#'
#' @param x A lisst object.
#'
#' @return A \code{data.frame}, with all columns of class \code{quantities} in 
#' exception of those that contain time information.
#'
#' @export

drop_lisst <- function(x) {
	structure(x, type = NULL, lproc = NULL, linst = NULL, 
			lmodl = NULL, zscat = NULL, class = "data.frame")
}


