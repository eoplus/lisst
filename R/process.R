# TODO
# add an subset and apply function that preserves structure, units and attributes.

lgdate <- function(lraw, yr) {

	if(!is(lraw, "lisst"))
		stop("lraw must be a lisst object", call. = FALSE)
	if(!(attr(lraw, "type") == "raw" || attr(lraw, "type") == "cor"))
		stop("lraw must be a lisst object of type raw or cor", call. = FALSE)

	if(attr(lraw, "linst")$mod == "100") {
		julian <- floor(lraw[, 39] / 100)

		if(missing(yr)) {
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
        	yr <- rep(yr, nrow(lraw))
		if(length(id) > 0)
			for(i in 1:length(id)) {
				yr[id[i]:length(yr)] <- yr[id[i]:length(yr)]+1
			}
		hour   <- round(((lraw[, 39] / 100) - julian) * 100)
		min    <- floor(lraw[, 40] / 100)
        	sec    <- round(((lraw[, 40] / 100) - min) * 100)
        	dates  <- as.POSIXct(paste(yr, julian, hour, min, sec, sep = "-"), 
			format = "%Y-%j-%H-%M-%S")
		return(dates)
	}

	if(attr(lraw, "linst")$mod == "200") {
		cat("To be added soon...")
	}
}



#' Retrieve the Particle Volume Scattering Function (VSF) from calibrated LISST data
#'
#' @param lcal object of class lisst and type cal.
#'
#' References:
#'
# Agrawal, Yogesh C. 2005. The optical volume scattering function: Temporal and vertical 
# variability in the water column off the New Jersey coast. Limnology and Oceanography 50, 6, 
# 1787-1794. DOI: 10.4319/lo.2005.50.6.1787
#

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
