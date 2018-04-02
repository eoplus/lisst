#' Subsenting of lisst objects
#' 
#' Performs subsetting of keeping the attributes of the lisst object.
#'
#' @export

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

#' Retrieve the LISST corrected digital counts
#'
#' The function retrieve the corrected digital counts (lisst object type 'cor')
#' from raw digital counts ('raw') or from calibrated values ('cal').
#'
#' @param lo A lisst object of type 'raw', 'cor' or 'cal'.
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
#' @export

lgetcor <- function(lo) {
	typ <- attr(lo, "type")
	if(typ == 'cor')
		return(lo)
	if(!(typ == 'cal' || typ == 'raw'))
		stop("Corrected digital counts can only be retrieved from a 'raw' or 'cal' lisst object", call. = FALSE)
	zscat <- attr(lo, "zscat")
	if(is.na(zscat[1]))
		stop("zscat data is missing from lisst object. Run read_lisst with zscat file path", call. = FALSE)
	if(is.na(attr(lo, "sn")))
		stop("Corrected digital counts requires instrument specific information. Run read_list with sn", call. = FALSE)

	linst <- attr(lo, "linst")
	lmodl <- attr(lo, "lmodl")
	if(typ == 'raw') {
		aw670 <- 0.439
		bw670 <- 0.0005808404
		wext <- exp(-(aw670 + bw670) * as.numeric(lmodl$pl))

		tau <- lo[, "Laser transmission"] * zscat["Laser reference"] / 
			zscat["Laser transmission"] / lo[, "Laser reference"]

		lo[, 1:lmodl$nring] <- lo[, 1:lmodl$nring] / tau - (rep(zscat[1:lmodl$nring], 
			each = nrow(lo)) * lo[, "Laser reference"] / zscat["Laser reference"])
		lo[, 1:lmodl$nring] <- lo[, 1:lmodl$nring] / wext
		lo[, 1:lmodl$nring] <- rep(linst$ringcf, each = nrow(lo)) * lo[, 1:lmodl$nring]
		lo[, 1:lmodl$nring][lo[, 1:lmodl$nring] < 0] <- 0
	} else {
		lo[, "Laser transmission"]   <- (lo[, "Laser transmission"] - linst$lpowcc[2]) / linst$lpowcc[1]
		lo[, "Battery voltage"]      <- (lo[, "Battery voltage"] - linst$battcc[2]) / linst$battcc[1]
		lo[, "External input 1"]     <- (lo[, "External input 1"] - linst$extrcc[2]) / linst$extrcc[1]
		lo[, "Laser reference"]      <- (lo[, "Laser reference"] - linst$lrefcc[2]) / linst$lrefcc[1]
		lo[, "Depth"]                <- (lo[, "Depth"] - linst$dpthcc[2]) / linst$dpthcc[1]
		lo[, "Temperature"]          <- (lo[, "Temperature"] - linst$tempcc[2]) / linst$tempcc[1]
		lo[, 1:lmodl$nring] <- as.data.frame(set_units(as.matrix(lo[, 1:lmodl$nring]), mW) / linst$ringcc)
		lo <- lo[, -which(names(lo) == "Optical transmission")]
		lo <- lo[, -which(names(lo) == "Beam attenuation")]
		lo <- lo[, -which(names(lo) == "Time")]
	}
	attr(lo, "type")  <- "cor"
	return(lo)
}

#' Retrieve the LISST calibrated values
#'
#' The function retrieve the calibrated values (lisst object type 'cal') from 
#' raw digital counts ('raw') or from corrected digital counts ('cor').
#'
#' @param lo A lisst object of type 'raw', 'cor' or 'cal'.
#'
#' @details When supplying a lisst 'raw' object, it must contain information on 
#' background values and instrument specific calibration constants. This is the
#' case when the 'raw' lisst object was created with zscat and sn arguments set
#' through \code{read_lisst}.
#'
#' Note that the values are already de-attenuated from pure water extinction in
#' 'cor' generation (precursor to 'cal').
#'
#' @export

lgetcal <- function(lo, yr) {
	typ <- attr(lo, "type")
	if(typ == 'cal')
		return(lo)
	if((!typ == 'raw' || typ == 'cor'))
		stop("Calibrated units can only be retrieved from a 'raw' or 'cor' lisst object", call. = FALSE)
	if(typ == 'raw') lo <- lgetcor(lo)
	zscat <- attr(lo, "zscat")
	linst <- attr(lo, "linst")
	lmodl <- attr(lo, "lmodl")

	tau <- lo[, "Laser transmission"] * zscat["Laser reference"] / 
		zscat["Laser transmission"] / lo[, "Laser reference"]

	lo[, "Laser transmission"]   <- lo[, "Laser transmission"] * linst$lpowcc[1] + linst$lpowcc[2]
	lo[, "Battery voltage"]      <- lo[, "Battery voltage"] * linst$battcc[1] + linst$battcc[2]
	lo[, "External input 1"]     <- lo[, "External input 1"] * linst$extrcc[1] + linst$extrcc[2]
	lo[, "Laser reference"]      <- lo[, "Laser reference"] * linst$lrefcc[1] + linst$lrefcc[2]
	lo[, "Depth"]                <- lo[, "Depth"] * linst$dpthcc[1] + linst$dpthcc[2]
	lo[, "Temperature"]          <- lo[, "Temperature"] * linst$tempcc[1] + linst$tempcc[2]
	lo[, "Optical transmission"] <- tau
	lo[, "Beam attenuation"]     <- -log(tau) / lmodl$pl
	lo[, 1:lmodl$nring] <- as.data.frame(set_units(as.matrix(lo[, 1:lmodl$nring]) * linst$ringcc, mW))

	if(missing(yr)) yr <- NULL
	lo$Time <- lgdate(lo, yr)
	attr(lo, "type")  <- "cal"
	return(lo)
}

#' Retrieve VSF from LISST data
#'
#' The function retrieves the absolute particle Volume Scattering Function 
#' (1/m/sr) for LISST data.
#'
#' @param lo A lisst object of type 'raw', 'cor' or 'cal'.
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
#' @export

lgetvsf <- function(lo) {

	if(!is(lo, "lisst"))
		stop("lo must be a lisst object", call. = FALSE)
	if(attr(lo, "type") != "cal")
		stop("lo must be a lisst object of type cal", call. = FALSE)

	linst <- attr(lo, "linst")
	lmodl <- attr(lo, "lmodl")

	if(linst$mod == "100") {
		wang  <- c(lmodl$wang[1, 2], lmodl$wang[, 1])
		zscat <- attr(lo, "zscat")

		for(i in 1:32)
			lo[, i] <- set_units(lo[, i] / lo[, 36] / (pi * lmodl$pl * 
				(wang[i]^2 - wang[i+1]^2) / 6), 1/m/sr)
		attr(lo, "type") <- "vsf"
		return(lo)
	}
	if(ancd$mod == "m200") {
		cat("To be added soon...")
	}
}

#' Retrieve PSD in number concentration
#'
#' The function converts the PSD in volume concentration (µL/L, ppm) to number 
#' concentration (particle/L/µm).
#'
#' @param lo A lisst object of type 'vol'.
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
#' @export

lgetpnc <- function(lo) {
	if(!is(lo, "lisst"))
		stop("lo must be a lisst object", call. = FALSE)
	typ <- attr(lo, "type")
	if(typ == "pnc")
		return(lo)
	if(typ != "vol")
		stop("lo must be a lisst object of type 'vol'", call. = FALSE)

	linst <- attr(lo, "linst")
	lmodl <- attr(lo, "lmodl")

        bins  <- lmodl$binr[[linst$ity]]
	nconc <- set_units(4 * pi * (bins[, 3] / 2)^3 / 3, L)
        binl  <- bins[, 2] - bins[, 1]
	fact  <- 1 / nconc / binl
	for(i in 1:lmodl$nring) {
		lo[, i] <- units::set_units(lo[, i] * fact[i] , 1/L/µm)
	}
	attr(lo, "type") <- "pnc"
	return(lo)
}

#' Retrieve PSD in volume concentration
#'
#' The function converts the PSD in particle concentration (particle/L/µm) back 
#' to the original data in volume concentration (µL/L, ppm).
#'
#' @param lo A lisst object of type 'pnc'.
#'
#' @details It merelly reverts the multiplication factors used by \code{lgetpnc}.
#'
#' @export

lgetvol <- function(lo) {
	if(!is(lo, "lisst"))
		stop("lo must be a lisst object", call. = FALSE)
	typ <- attr(lo, "type")
	if(typ == "vol")
		return(lo)
	if(typ != "pnc")
		stop("lo must be a lisst object of type 'pnc'", call. = FALSE)

	linst <- attr(lo, "linst")
	lmodl <- attr(lo, "lmodl")

        bins  <- lmodl$binr[[linst$ity]]
	nconc <- set_units(4 * pi * (bins[, 3] / 2)^3 / 3, L)
        binl  <- bins[, 2] - bins[, 1]
	fact  <- 1 / nconc / binl
	for(i in 1:lmodl$nring) {
		lo[, i] <- units::set_units(lo[, i] / fact[i] , ppm)
	}
	attr(lo, "type") <- "vol"
	return(lo)
}

#' Get LISST measurement dates
#' 
#' The function retrieves the date/time of measurements for the records in the 
#' lisst object as POSIXct objects. The time zone of the dates will depend on 
#' user tz options.
#'
#' @param lo A lisst object of type 'raw', 'cor' or 'cal'.
#' @param yr The year of first measurement in the lisst object. Can be omitted. 
#'           Ignored for LISST-200X.
#'
#' @details If yr is not provided (or is NULL) teh function will perform a 'best 
#' guess' (with a warning) based on the curent date and the julian date of first 
#' measurement. Simply, if the julian day of first measurement is greater than
#' the curent julian day, the previous from current year is assumed and the 
#' current year otherwise.
#'
#' @export

lgdate <- function(lo, yr) {
	if(!is(lo, "lisst"))
		stop("lo must be a lisst object", call. = FALSE)
	typ <- attr(lo, "type")
	if((typ == "cal" || typ == "vol" || typ == "pnc") && (sum(names(lo) == "Time") > 0))
		return(lo[, "Time"])

	if(attr(lo, "lmodl")$mod == "100") {
		julian <- floor(lo[, 39] / 100)
		if(missing(yr) || is.null(yr)) {
			ct <- Sys.time()
			cj <- as.numeric(format(ct, "%j"))
			yr <- as.numeric(format(ct, "%Y"))
			if(julian[1] > cj) {
				yr <- yr - 1
				warning("yr not provided; assuming previous year based on",
					" julian dates", call. = FALSE)
			} else {
				warning("yr not provided, assuming current year based on", 
					" julian dates", call. = FALSE)
			}
		}
		id <- which(diff(julian) <= -364)
        	yr <- rep(yr, nrow(lo))
		if(length(id) > 0)
			for(i in 1:length(id)) {
				yr[id[i]:length(yr)] <- yr[id[i]:length(yr)]+1
			}
		hour   <- round(((lo[, 39] / 100) - julian) * 100)
		min    <- floor(lo[, 40] / 100)
        	sec    <- round(((lo[, 40] / 100) - min) * 100)
        	dates  <- as.POSIXct(paste(yr, julian, hour, min, sec, sep = "-"), 
			format = "%Y-%j-%H-%M-%S")
		return(dates)
	}

	if(attr(lo, "linst")$mod == "200") {
		cat("To be added soon...")
	}
}

