# TODO
# add an subset and apply function that preserves structure, units and attributes.

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
		tau <- lo[, "Laser transmission"] * zscat["Laser reference"] / 
			zscat["Laser transmission"] / lo[, "Laser reference"]
		lo[, 1:lmodl$nring] <- lo[, 1:lmodl$nring] / tau - (rep(zscat[1:lmodl$nring], 
			each = nrow(lo)) * lo[, "Laser reference"] / zscat["Laser reference"])
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
					"julian dates", call. = FALSE)
			} else {
				warning("yr not provided, assuming current year based on", 
					"julian dates", call. = FALSE)
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



#' Retrieve the Particle Volume Scattering Function (VSF) from calibrated LISST 
#' data.
#'
#' @param lcal object of class lisst and type cal.
#'
#' @references
#'
#' Agrawal, Yogesh C. 2005. The optical volume scattering function: Temporal and 
#' vertical variability in the water column off the New Jersey coast. Limnology 
#' and Oceanography 50, 6, 1787-1794. DOI: 10.4319/lo.2005.50.6.1787
#'
#' @export

getvsf <- function(lcal) {

	if(!is(lcal, "lisst"))
		stop("lcal must be a lisst object", call. = FALSE)
	if(attr(lcal, "type") != "cal")
		stop("lcal must be a lisst object of type cal", call. = FALSE)

	linst <- attr(lcal, "linst")
	lmodl <- attr(lcal, "lmodl")

	if(linst$mod == "100") {
		wang  <- c(lmodl$wang[1, 2], lmodl$wang[, 1])
		zscat <- attr(lcal, "zscat")
		aw670 <- 0.439
		bw670 <- 0.0005808404
		wext <- exp(-(aw670 + bw670) * as.numeric(lmodl$pl))

		for(i in 1:32)
			lcal[, i] <- set_units(lcal[, i] / lcal[, 36] / (wext * pi * lmodl$pl * 
				(wang[i]^2 - wang[i+1]^2) / 6), 1/m/sr)
		attr(lcal, "type") <- "vsf"
		return(lcal)
	}
	if(ancd$mod == "m200") {
		cat("To be added soon...")
	}
}

#
# Function: getnconc
#
