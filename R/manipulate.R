
# Get lisst specific attributes
# This was mainly created for use in the c.lisst function.

.lattributes <- function(x) {
	attributes(x)[c('type', 'lproc', 'linst', 'lmodl', 'zscat')]
}

#' Concatenate lisst objects
#' 
#' This function allows to concatenate two or more lisst objects into a single
#' object. Checks are made to ensure the final object is coherent (same 
#' instrument and inversion type, if applicable) and conversion between types is
#' made as necessary. See details.
#'
#' @param ... lisst objects to concatenate.
#' 
#' @details The first lisst object is taken as reference, i.e., all subsequent 
#' lisst objects will be converted to its type ('raw', 'cor', 'cal', 'vsf', 
#' 'vol' or 'pnc') before the concatenation. A check on the lisst attributes is 
#' made to ensure consistency that all data is from a single instrument/model 
#' and from the same inversion type (if applicable).
#'
#' @return A lisst object with type equal to the first argument and containing 
#' data from all arguments.
#'
#' @examples
#' l1 <- donkmeer_bin
#' l2 <- lget(donkmeer_bin, 'raw')
#' l3 <- c(l1, l2)
#' 
#' @export

c.lisst <- function(...) {
	x   <- list(...)
	for(i in 1:length(x)) {
		if(!is.lisst(x[[1]])) stop(paste("Argument", i, "is not a lisst object"), call. = F)
	}
	typ <- attr(x[[1]], "type")
	x   <- lapply(x, lget, type = typ)
	for(i in 2:length(x)) {
		if(!identical(.lattributes(x[[1]]), .lattributes(x[[i]])))
			stop("Concatenation of lisst objects requires them to be from the same ",
				"instrument/model and the same inversion type (if applicable)", 
				call. = FALSE)
	}
	do.call(rbind, x)
}

#' Subset lisst objects
#' 
#' Performs subsetting by index, depth or time.
#'
#' @param i An integer vector, a time object (e.g., POSIXct), or a character 
#' string specifying deph range (in meters) or time range in the ISO 8601 
#' format.
#' @param j An integer vector or a character string to macth a column name.
#'
#' @details Subseting can be done by sample index using a vector of integers, 
#' by depth using a string of the form 'depth1|depth2', or by time using a time 
#' object vector (e.g. POSIXct) or a caracter string in the format of ISO 8601.
#' The time subsetting actully uses the xts subsetting functions, see 
#' \code{?`[.xts`} for details.
#'
#' @examples
#' donkmeer_bin[1:3, ]         # First three samples
#' donkmeer_bin['0|5', ]       # First five meters
#' donkmeer_bin['2018-03', ]   # Samples for march 2018
#' donkmeer_bin['2017-06/2018-06', ] # Samples between June 2017 to June 2018
#'
#' @seealso \code{`[.xts`}
#'
#' @export

`[.lisst` <- function(x, i, j, ..., drop = TRUE) {
	stopifnot(is.lisst(x))
	lxts <- xts(rep(1, nrow(x)), order.by = x$Time, unique = FALSE)
	if(!missing(i)) {
		if(is.character(i)) {
			if(length(i) > 1) {
			# Will be used when station names are recorded. Now will just throw an error.
			stop("Character indexing for lisst objects must be length one", call. = FALSE)
			} else if(length(grep("|", i, fixed = TRUE)) > 0) {
				range <- as.numeric(unlist(strsplit(i, "|", fixed = TRUE)))
				depth <- drop_quantities(x$Depth)
				i <- which(depth >= range[1] & depth < range[2])
			} else {
				i <- xts:::`[.xts`(lxts, i, which.i = TRUE)
			}
		} else if(xts::timeBased(i)) {
			i <- xts:::`[.xts`(lxts, i, which.i = TRUE)
		}
	}
	if(((missing(i) && length(j) == 1) || (missing(j) && length(i) == 1)) && drop) {
		#x <- drop_lisst(x)
		x <- NextMethod(drop = drop)
	} else {
		x <- NextMethod(drop = drop)
		lxts <- lxts[i, ]
		structure(x, 
			"type"  = attr(x, "type"),
			"lproc" = attr(x, "lproc"),
			"linst" = attr(x, "linst"),
			"lmodl" = attr(x, "lmodl"),
			"zscat" = attr(x, "zscat"),
#			"lxts"  = lxts,
			class   = c("lisst", "data.frame"))
	}
}

#' Convert between LISST data representations
#'
#' This set of functions allow for the convertion between lisst data types.
#'
#' @param x    A lisst object.
#' @param type A character specifing the lisst data type. One of 'raw', 'cor', 
#' 'cal', 'vsf', 'vol', 'pnc'. See details.
#' 
#' @details lget dispatch the appropriate lgetxxx function acording to parameter 
#' type. As of this version is only possible to convert between 'raw', 'cor', 
#' 'cal' and 'vsf' or between 'vol' and 'pnc' (inversion not implemented). 
#' \describe{
#'   \item{raw}{The raw digital counts as recorded by the LISST instrument.}
#'   \item{cor}{The corrected digital counts, i.e., the raw counts de-attenuated 
#'              for the particle \strong{and water} extinction, background 
#'              subtracted and compensated for area deviations from nominal 
#'              values.}
#'   \item{cal}{The calibrated values, i.e., the instrument specific calibration 
#'              constants applied to the corrected values (for all variables). 
#'              Aditionally, the transmittance due to particles and the particle 
#'              beam attenuation are added to the lisst object.}
#'   \item{vsf}{The volume scattering function, i.e., the calibrated values 
#'              normalized to incident power, the solid angle of the detectors 
#'              and the path of water generating the signal.}
#'   \item{vol}{The particle volume concentration (ppm volume) per size bin, as 
#'		inverted from the scattering data by the LISST-SOP.}
#'   \item{pnc}{The particle number concentration per volume and µm size.}
#' }
#' 
#' See the documentation of the lget functions for further details.
#'
#' @examples
#' l_vsf <- donkmeer_bin
#' l_cal <- lget(l_vsf, 'cal')
#' l_cor <- lget(l_vsf, 'cor')
#' l_raw <- lget(l_vsf, 'raw')
#' l_vsf <- lget(l_raw, 'vsf')
#'
#' @seealso \code{\link{lgetraw}}, \code{\link{lgetcor}}, \code{\link{lgetcal}},
#' \code{\link{lgetvsf}}, \code{\link{lgetvol}}, \code{\link{lgetpnc}}
#'
#' @export

lget <- function(x, type) {
	stopifnot(is.lisst(x))
	switch(type, 
		"raw" = lgetraw(x),
		"cor" = lgetcor(x),
		"cal" = lgetcal(x),
		"vsf" = lgetvsf(x),
		"vol" = lgetvol(x),
		"pnc" = lgetpnc(x),
		stop("type must be one of: 'raw', 'cor', 'cal', 'vsf', 'vol', 'pnc'", call. = FALSE)
	)
}

#' Retrieve the LISST raw digital counts
#'
#' The function retrieve the raw digital counts (lisst object type 'raw')
#' from corrected digital counts ('cor'), calibrated values ('cal'), or  
#' volume scattering function ('vsf').
#'
#' @param x A lisst object of type 'raw', 'cor', 'cal' or 'vsf'.
#'
#' @examples
#' l_raw <- lgetraw(donkmeer_bin)
#'
#' @seealso \code{\link{lget}} \code{\link{lgetcor}}, \code{\link{lgetval}}, 
#' \code{\link{lgetvsf}}, \code{\link{lgetvol}}, \code{\link{lgetpnc}}
#'
#' @export

lgetraw <- function(x) {
	stopifnot(is.lisst(x))
	typ <- attr(x, "type")
	if(typ == 'raw')
		return(x)
	else if(typ == 'cal' || typ == 'vsf') x <- lget(x, 'cor')
	else if(typ != 'cor')
		stop("Raw digital counts can only be retrieved from a 'cor', 'cal' or 'vsf' lisst ",
			"object", call. = FALSE)

	zscat <- attr(x, "zscat")
	linst <- attr(x, "linst")
	lmodl <- attr(x, "lmodl")

	wext <- quantities::set_quantities(exp(-(aw670 + bw670) * as.numeric(lmodl$pl)), 1, 0)
	tau  <- x[, "Laser transmission"] * zscat[, "Laser reference"] / 
		zscat[, "Laser transmission"] / x[, "Laser reference"]

	for(i in 1:lmodl$nring) {
		x[, i] <- x[, i] / linst$ringcf[i]
		x[, i] <- x[, i] * wext
		x[, i] <- (x[, i] + (zscat[, i] * x[, "Laser reference"] / 
			zscat[, "Laser reference"])) * tau
	}

	attr(x, "type")  <- "raw"
	x
}

#' Retrieve the LISST corrected digital counts
#'
#' The function retrieve the corrected digital counts (lisst object type 'cor')
#' from raw digital counts ('raw'), calibrated values ('cal') or volume 
#' scattering function ('vsf').
#'
#' @param x A lisst object of type 'raw', 'cor', 'cal' or 'vsf'.
#'
#' @details When supplying a lisst 'raw' object, it must contain information on 
#' background values and instrument specific calibration constants. This is the
#' case when the 'raw' lisst object was created with zscat and sn arguments set
#' through \code{read_lisst}.
#'
#' For consistency with LISST processed data provided by manufacturer, the 
#' optical transmission and the beam attenuation do not include the effect of 
#' pure water. However, the digital counts in the ring detectors are corrected 
#' for the aditional transmission loss due to pure water absorption and 
#' scattering. Values for the pure water absorption and scattering at 670 nm are 
#' taken from the Water Optical Properties Processor (WOPP) and correspond to 
#' 0.439 and 5.808e-4 1/m, respectivelly.
#'
#' @seealso \code{\link{lget}} \code{\link{lgetraw}}, \code{\link{lgetval}}, 
#' \code{\link{lgetvsf}}, \code{\link{lgetvol}}, \code{\link{lgetpnc}}
#'
#' @examples
#' l_cor <- lgetcor(donkmeer_bin)
#'
#' @export

lgetcor <- function(x) {
	stopifnot(is.lisst(x))
	typ <- attr(x, "type")
	if(typ == 'cor') return(x)
	else if(typ == 'vsf') {
		x <- lget(x, 'cal')
		typ <- attr(x, "type")
	}
	else if(!(typ == 'cal' || typ == 'raw'))
		stop("Corrected digital counts can only be retrieved from a 'raw', 'cal' or 'vol' ",
			"lisst object", call. = FALSE)
	zscat <- attr(x, "zscat")
	if(is.na(zscat[1]))
		stop("zscat data is missing from lisst object. Run read_lisst with zscat file ",
			"path set", call. = FALSE)
	if(is.na(attr(x, "linst")$sn))
		stop("Corrected digital counts requires instrument specific information. Run ",
			"read_list with sn set", call. = FALSE)

	linst <- attr(x, "linst")
	lmodl <- attr(x, "lmodl")
	if(typ == 'raw') {
		wext <- set_quantities(exp(-(aw670 + bw670) * as.numeric(lmodl$pl)), 1, 0)
		tau  <- x[, "Laser transmission"] * zscat[, "Laser reference"] / 
			zscat[, "Laser transmission"] / x[, "Laser reference"]

		for(i in 1:lmodl$nring) {
			x[, i] <- x[, i] / tau - (zscat[, i] * x[, "Laser reference"] / 
				zscat[, "Laser reference"])
			x[, i] <- x[, i] / wext
			x[, i] <- x[, i] * linst$ringcf[i]
			x[, i][drop_quantities(x[, i]) < 0] <- .zero
		}
	} else {
		x[, "Laser transmission"] <- (x[, "Laser transmission"] - linst$lpowcc[2]) / linst$lpowcc[1]
		x[, "Battery voltage"]    <- (x[, "Battery voltage"]    - linst$battcc[2]) / linst$battcc[1]
		x[, "External input 1"]   <- (x[, "External input 1"]   - linst$extrcc[2]) / linst$extrcc[1]
		x[, "Laser reference"]    <- (x[, "Laser reference"]    - linst$lrefcc[2]) / linst$lrefcc[1]
		x[, "Depth"]              <- (x[, "Depth"]              - linst$dpthcc[2]) / linst$dpthcc[1]
		x[, "Temperature"]        <- (x[, "Temperature"]        - linst$tempcc[2]) / linst$tempcc[1]
		for(i in 1:lmodl$nring) x[, i] <- x[, i] / linst$ringcc
		x <- x[, -which(names(x) == "Optical transmission"), drop = FALSE]
		x <- x[, -which(names(x) == "Beam attenuation"), drop = FALSE]
	}
	attr(x, "type")  <- "cor"
	x
}

#' Retrieve the LISST calibrated values
#'
#' The function retrieve the calibrated values (lisst object type 'cal') from 
#' raw digital counts ('raw'), corrected digital counts ('cor') or volume 
#' scattering function ('vsf').
#'
#' @param x A lisst object of type 'raw', 'cor', 'cal' or 'vsf'.
#'
#' @details When supplying a lisst 'raw' object, it must contain information on 
#' background values and instrument specific calibration constants. This is the
#' case when the 'raw' lisst object was created with zscat and sn arguments set
#' through \code{read_lisst}.
#'
#' @seealso \code{\link{lget}} \code{\link{lgetraw}}, \code{\link{lgetcor}}, 
#' \code{\link{lgetvsf}}, \code{\link{lgetvol}}, \code{\link{lgetpnc}}
#'
#' @examples
#' l_cal <- lgetcal(donkmeer_bin)
#'
#' @export

lgetcal <- function(x) {
	stopifnot(is.lisst(x))
	typ <- attr(x, "type")
	if(typ == 'cal') return(x)
	else if(typ == 'raw') {
		x <- lget(x, 'cor')
		typ <- attr(x, "type")
	}
	else if(!(typ == 'cor' || typ == 'vsf'))
		stop("Calibrated units can only be retrieved from a 'raw', 'cor' or 'vsf' lisst ",
			"object", call. = FALSE)

	zscat <- attr(x, "zscat")
	linst <- attr(x, "linst")
	lmodl <- attr(x, "lmodl")

	if(typ == 'cor') {
		tau <- x[, "Laser transmission"] * zscat[, "Laser reference"] / 
			zscat[, "Laser transmission"] / x[, "Laser reference"]

		x[, "Laser transmission"]   <- x[, "Laser transmission"] * linst$lpowcc[1] + linst$lpowcc[2]
		x[, "Battery voltage"]      <- x[, "Battery voltage"]    * linst$battcc[1] + linst$battcc[2]
		x[, "External input 1"]     <- x[, "External input 1"]   * linst$extrcc[1] + linst$extrcc[2]
		x[, "Laser reference"]      <- x[, "Laser reference"]    * linst$lrefcc[1] + linst$lrefcc[2]
		x[, "Depth"]                <- x[, "Depth"]              * linst$dpthcc[1] + linst$dpthcc[2]
		x[, "Temperature"]          <- x[, "Temperature"]        * linst$tempcc[1] + linst$tempcc[2]
		x[, "Optical transmission"] <- tau
		x[, "Beam attenuation"]     <- set_units(drop_units(-log(tau) / lmodl$pl), 1/m)
		for(i in 1:lmodl$nring) x[, i] <- set_units(x[, i] * linst$ringcc, µW)
	} else {
		wang  <- as.errors(c(lmodl$wang[1, 2], lmodl$wang[, 1]))
		for(i in 1:lmodl$nring) 
			x[, i] <- set_units(x[, i] * x[, "Laser reference"] * (set_quantities(pi, 1, 0) * lmodl$pl * 
				set_units(wang[i]^2 - wang[i+1]^2, sr) / set_quantities(6, 1, 0)), µW)
	}

	attr(x, "type")  <- "cal"
	x
}

#' Retrieve VSF from LISST data
#'
#' The function retrieves the particle Volume Scattering Function (1/m/sr) for 
#' LISST data.
#'
#' @param x A lisst object of type 'raw', 'cor', 'cal' or 'vsf'.
#'
#' @details Types 'raw' and 'cor' are converted to 'cal' first. The function 
#' then normalizes the power (mW) measured by the ring detectors by their solid 
#' angle (sr), the path of water generating the signal (m), and the power 
#' entering the path (mW). Since when generating 'cor' the measured signal is
#' already de-attenuated from pure water extinction, no additional correction is 
#' necessary.
#'
#' @references
#' Agrawal, Yogesh C. 2005. The optical volume scattering function: Temporal and 
#' vertical variability in the water column off the New Jersey coast. Limnology 
#' and Oceanography 50, 6, 1787-1794. DOI: 10.4319/lo.2005.50.6.1787
#'
#' @seealso \code{\link{lget}} \code{\link{lgetraw}}, \code{\link{lgetcor}}, 
#' \code{\link{lgetcal}}, \code{\link{lgetvol}}, \code{\link{lgetpnc}}
#'
#' @examples
#' l_vsf <- lget(lgetraw(donkmeer_bin), 'vsf')
#' all.equal(l_vsf, donkmeer_bin)
#'
#' @export

lgetvsf <- function(x) {
	stopifnot(is.lisst(x))
	typ <- attr(x, "type")
	if(typ == 'vsf') return(x)
	else if(typ == 'raw' || typ == 'cor') {
		x <- lget(x, 'cal')
		typ <- attr(x, "type")
	}
	else if(typ != "cal")
		stop("VSF can only be retrieved from a 'raw', 'cor' or 'cal' lisst object", 
			call. = FALSE)

	linst <- attr(x, "linst")
	lmodl <- attr(x, "lmodl")
	zscat <- attr(x, "zscat")

	wang  <- as.errors(c(lmodl$wang[1, 2], lmodl$wang[, 1]))

	for(i in 1:lmodl$nring)
		x[, i] <- set_units(x[, i] / x[, "Laser reference"] / (set_quantities(pi, 1, 0) * lmodl$pl * 
				(wang[i]^2 - wang[i+1]^2) / set_quantities(6, 1, 0)), 1/m/sr)
	attr(x, "type") <- "vsf"
	x
}

#' Retrieve PSD in number concentration
#'
#' The function converts the PSD in volume concentration (µL/L, ppm) to number 
#' concentration (particle/L/µm).
#'
#' @param x A lisst object of type 'vol'.
#'
#' @details Volume concentration is converted to number concentration by using 
#' the volume of a sphere with radius equal to half the median particle size
#' for each bin. The number concentration is then the sphere equivalent number
#' concentration. The absolute magnitute of the PSD will be only approximate if 
#' particles are not spherical, but as long as the particles are not expected to
#' significantly change shape with size, the slope of the distribution will be 
#' accurate.
#'
#' @references
#' Buonassissi, C. J. and Dierssen, H. M. 2010. A regional comparison of 
#' particle size distributions and the power law approximation in oceanic and 
#' estuarine surface waters. J. Geophys. Res., 115, C10028. 
#' DOI:10.1029/2010JC006256.
#'
#' @seealso \code{\link{lget}} \code{\link{lgetraw}}, \code{\link{lgetcor}}, 
#' \code{\link{lgetcal}}, \code{\link{lgetvsf}}, \code{\link{lgetpnc}}
#'
#' @examples
#' l_pnc <- lgetpnc(donkmeer_pro)
#'
#' @export

lgetpnc <- function(x) {
	stopifnot(is.lisst(x))
	typ <- attr(x, "type")
	if(typ == "pnc") return(x)
	if(typ != "vol")
		stop("Particle number concentration can only be retrieved a 'vol' lisst object", 
			call. = FALSE)

	linst <- attr(x, "linst")
	lmodl <- attr(x, "lmodl")
	lproc <- attr(x, "lproc")

        bins  <- lmodl$binr[[lproc$ity]]
	nconc <- set_quantities(4 * pi * (bins[, 3] / 2)^3 / 3, L, 0)
        binl  <- set_errors(bins[, 2] - bins[, 1], 0)
	fact  <- nconc^-1 * binl^-1
	for(i in 1:lmodl$nring) {
		x[, i] <- set_units(x[, i] * fact[i] , 1/L/µm)
	}
	attr(x, "type") <- "pnc"
	x
}

#' Retrieve PSD in volume concentration
#'
#' The function converts the PSD in particle concentration (particle/L/µm) back 
#' to the original data in volume concentration (µL/L, ppm).
#'
#' @param x A lisst object of type 'pnc'.
#'
#' @details It merelly reverts the multiplication factors used by \code{lgetpnc}.
#'
#' @seealso \code{\link{lget}} \code{\link{lgetraw}}, \code{\link{lgetcor}}, 
#' \code{\link{lgetcal}}, \code{\link{lgetvsf}}, \code{\link{lgetvol}}
#'
#' @examples
#' l_vol <- lget(lgetpnc(donkmeer_pro), 'vol')
#' all.equal(l_vol, donkmeer_pro)
#'
#' @export

lgetvol <- function(x) {
	stopifnot(is.lisst(x))
	typ <- attr(x, "type")
	if(typ == "vol") return(x)
	if(typ != "pnc")
		stop("Particle volume concentration can only be retrieved from a 'pnc' lisst ",
			"object", call. = FALSE)

	linst <- attr(x, "linst")
	lmodl <- attr(x, "lmodl")
	lproc <- attr(x, "lproc")

        bins  <- lmodl$binr[[lproc$ity]]
	nconc <- set_quantities(4 * pi * (bins[, 3] / 2)^3 / 3, L, 0)
        binl  <- set_errors(bins[, 2] - bins[, 1], 0)
	fact  <- nconc^-1 * binl^-1
	for(i in 1:lmodl$nring) {
#		x[, i] <- set_units(x[, i] / fact[i], ppm) # Error; seems a bug in the development version of units...
		x[, i] <- set_units(drop_units(x[, i] / fact[i]) * 1e6, ppm)
	}
	attr(x, "type") <- "vol"
	return(x)
}

