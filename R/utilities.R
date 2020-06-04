
#' Extract metadata from lisst objects
#' 
#' This collection of functions allows to extract specific metadata stored in the lisst
#' object. 
#'
#' @details These functions are the recommended way to access the object's metadata since 
#' metadata organization might be subject to change. 
#'
#' @return The returned data is dependent on the specific metadata requested:
#'
#' \itemize{
#'   \item ltyi   - A character indicating the LISST SOP inversion type (rs or ss) if any;
#'   \item sizes  - A numeric matrix with the minimum, maximum and median sizes in each bin;
#'   \item angles - A numeric matrix with the minimum, maximum and median angles in each bin;
#'   \item ltime  - A POSIXct vector of recording times;
#'   \item zscat  - A numeric vector with background scattering values (blank);
#'   \item fzscat - A numeric vector with factory background scattering values (blank).
#' }
#'
#' @name metadata
NULL

#' @describeIn metadata Extract background scattering
#'
#' @param x A lisst object.
#'
#' @export

zscat  <- function(x) {
	stopifnot(is.lisst(x))
	return(.lattributes(x)$zscat)
}

#' @describeIn metadata Extract factory background scattering
#'
#' @export

fzscat <- function(x) {
	stopifnot(is.lisst(x))
	return(.LISSTi[[.lattributes(x)$linst$sn]])
}

#' @describeIn metadata Extract inversion type
#'
#' @export

lity   <- function(x) {
	stopifnot(is.lisst(x))
	return(.lattributes(x)$lproc$ity)
}

#' @describeIn metadata Extract lisst time index
#'
#' @export

ltime <- function(x) {
	stopifnot(is.lisst(x))
        rn <- rownames(x)
	ms <- min(nchar(rn)[nchar(rn) > 20])
        id <- which(nchar(rn) < 20)
        rn[id] <- NA
        rn <- substring(rn, 1, ms)
	as.POSIXct(rn)
}

#' @describeIn metadata Extract reference measuring angles
#'
#' @param units  Desired units.
#' @param medium Desired medium (water or air).
#'
#' @export

angles <- function(x, units, medium = c('water', 'air')) {
	stopifnot(is.lisst(x))
	angles <- .lattributes(x)$lmodl$wang
	if(medium[1] == 'air')
		asin(sin(angles) / 1.33)
	if(!missing(units)) units(angles) <- units
	return(angles)
}

#' @describeIn metadata Extract reference inversion sizes
#'
#' @export

sizes <- function(x, units) {
	stopifnot(is.lisst(x))
	sizes <- .lattributes(x)$lmodl$binr[[lity(x)]]
	if(!missing(units)) units(sizes) <- units
	return(sizes)
}


