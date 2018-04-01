

#
# Function: lisst_bin
#
# Read a LISST binary file for supported instruments. The function perform minimal tests of 
# input arguments and is not intended to be used directly by the user, but called by 
# read_lisst function.
#
# Input:
# fl    - path to binary file
# zscat - path to ASCII background data file
# pl    - path length in m. If not provided the system will assumed standard path length for
#         the instrument model, i.e., no path reduction module in use.
# sn    - serial number of the instrument
#
# Output:
# A lisstBin object where the number of columns depend on instrument type. The binary data is
# converted to calibrated units for all paramaters except scattering into ring detectors. For 
# the ring detectors the value returned is the "net scattering" in digital counts, i.e., 
# the scattering de-attenuated, background subtracted and corrected for area of the rings. 
# Those are the values required for input to the vendor proprietary inversion code to retrieve
# particle size distribution. The output is essentially equal to a processed ASCII file, since 
# the optical transmission and beam attenuation are also calculated, with the difference being 
# the variables of the ring detectors, being scattering for this binary read and PSD for an 
# ASCII file.
#
# Details:
# The return of this function is not expected to be directly usefull, but is provided as an
# intermediary step for further inversion to retrieve PSD or normalization to retrieve VSF.
#
# qck the quality contrl is ignored if out = "raw"

lisst_bin <- function(fl, zscat, sn, pl, yr, out = c("raw", "cor", "cal"), qck = TRUE) {
	sn    <- as.character(sn)
	linst <- .LISSTi[[sn]]
	lmodl <- .getmodp(linst)

	if(missing(pl))
		warning(paste("pl not provided; assuming standard path length of instrument"), 
			call. = FALSE)
	else 
			lmodl$pl <- set_units(pl, m)

	if(linst$mod == "100") {
		lraw <- readBin(fl, "integer", n = file.info(fl)$size, size = 1, signed = FALSE, 
			endian = "little")
		lraw <- lraw[seq(1, length(lraw), 2)] * 256 + lraw[seq(2, length(lraw), 2)]
		nrows <- floor(length(lraw) / lmodl$bnvar)
		lraw <- lraw[1:(nrows * lmodl$bnvar)]
		lraw <- matrix(lraw, nrow = nrows, byrow = TRUE)
		if(linst$X) {
			lraw[, 1:32] <- lraw[, 1:32] / 10
		}

		zscat <- as.numeric(read.table(zscat)[, 1])

		lraw <- as.data.frame(lraw)
		colnames(lraw) <- lmodl$lvarn[1:lmodl$bnvar]
		attr(lraw, "sn")    <- sn
		attr(lraw, "type")  <- "raw"
		attr(lraw, "linst") <- linst
		attr(lraw, "lmodl") <- lmodl
		attr(lraw, "zscat") <- zscat
		attr(lraw, "class") <- c("lisst", "data.frame")

		if(out == "raw")
			return(lraw)

		tau <- lraw[, 33] * zscat[36] / zscat[33] / lraw[, 36]
		lraw[, 1:32] <- lraw[, 1:32] / tau - (rep(zscat[1:32], each = nrow(lraw)) * lraw[, 36] / zscat[36])
		lraw[, 1:32] <- rep(linst$ringcf, each = nrow(lraw)) * lraw[, 1:32]
		lraw[, 1:32][lraw[, 1:32] < 0] <- 0

		attr(lraw, "type")  <- "cor"

		if(out == "cor")
			return(lraw)

		lraw[, 33] <- lraw[, 33] * linst$lpowcc[1] + linst$lpowcc[2]
		lraw[, 34] <- lraw[, 34] * linst$battcc[1] + linst$battcc[2]
		lraw[, 35] <- lraw[, 35] * linst$extrcc[1] + linst$extrcc[2]
		lraw[, 36] <- lraw[, 36] * linst$lrefcc[1] + linst$lrefcc[2]
		lraw[, 37] <- lraw[, 37] * linst$dpthcc[1] + linst$dpthcc[2]
		lraw[, 38] <- lraw[, 38] * linst$tempcc[1] + linst$tempcc[2]
		lraw[lmodl$lvarn[41]] <- tau
		lraw[lmodl$lvarn[42]] <- -log(tau) / lmodl$pl
		lraw[, 1:32] <- as.data.frame(set_units(as.matrix(lraw[, 1:32]) * linst$ringcc, mW))
		if(missing(yr))
			lraw[, "Time"] <- lgdate(lraw)
		else
			lraw[, "Time"] <- lgdate(lraw, yr)
		lraw <- lraw[, c(1:38, 41:43)]
		attr(lraw, "type")  <- "cal"
		return(lraw)
	}
	if(ancd$mod == "200") {
		cat("To be added soon...")
	}
}



#
# Function: lisst_asc
#
# Read a LISST processed ASCII file for supported instruments. The function perform minimal 
# tests of input arguments and is not intended to be used directly by the user, but called by 
# read_lisst function.
#
# Input:
# fl    - path to binary file
# pl    - path length in m. If not provided the system will assumed standard path length for
#         the instrument model, i.e., no path reduction module in use.
# sn    - serial number of the instrument
#
# Output:
# A lisst object where the number of columns depend on instrument type. Note that the units 
# of the volume concentration are ppm (parts per million by volume), which is equivalent to
# µL/L.
#
# Details:
# The return of this function is not expected to be directly usefull, but is provided as an
# intermediary step for further inversion to retrieve PSD or normalization to retrieve VSF.
#

lisst_asc <- function(fl, sn, pl = NULL) {
	sn <- paste0("sn_", sn)
	ancd  <- lisst_inst[[sn]]
	model <- lisst_model[[ancd$mod]]
	ldat <- read.table(fl, header = FALSE)
	colnames(ldat) <- svarnames[model$vars]
	if(is.null(pl)) {
		warning(paste("Beam attenuation calculated assuming standard path length",
			"of instrument"), call. = FALSE)
		mf <- 1
	} else {
		mf <-  model$pl / pl
	}
	ldat$Atten <- ldat$Atten * mf
	ldat <- as.data.frame(ldat)

	if(any(ldat$OTrans > 0 && ldat$OTrans < 0.1)) {
		warning(paste("Values for samples with optical transmittance lower than 0.1",
			"are set to zero"), call. = FALSE)
		id <- which(lraw$OTrans < 0.1)
		lraw[id, 1:model$nring] <- 0
		lraw$OTrans[id] <- 0
		lraw$Atten[id] <- 0
	}
	if(any(ldat$OTrans >= 0.1 && ldat$OTrans < 0.3))
		warning(paste("Optical transmittances below 0.3 might suffer from significant",
			"multiple scattering. Values bettewn 0.3 and 0.1 were kept but", 
			"caution is avised."), call. = FALSE)

	for(i in 1:model$nring)
		ldat[, i] <- ldat[, i] * varunits[[1]]
	for(i in (model$nring + 1):length(model$vars))
		units(ldat[, i]) <- varunits[[model$vars[i]-35]]
	colnames(ldat) <- lvarnames[model$vars]
	attr(ldat, "sn") <- sn
	ldat
}

read_lisst <- function(fl, sn, mode = c("binary", "ascii")) {

	if(is.null(lisst_inst[[paste("sn_", sn)]]))
		stop("Instrument not registered", call. = FALSE)
	if(!file.exists(fl))
		stop(paste("File", fl, "not found"), call. = FALSE)
	if(!file.exists(zscat))
		stop(paste("File", zscat, "not found"), call. = FALSE)



		if(qck && any(lraw$OTrans < 0.1)) {
			warning(paste("qck set to TRUE. Values for samples with optical", 
				"transmittance lower than 0.1 were set to zero"), call. = FALSE)
			id <- which(lraw$OTrans < 0.1)
			lraw[id, 1:model$nring] <- 0
			lraw$OTrans[id] <- 0
			lraw$Atten[id] <- 0
		}
		if(any(lraw$OTrans >= 0.1 && lraw$OTrans < 0.3))
			warning(paste("Optical transmittances below 0.3 might suffer from",
				"significant multiple scattering. Caution is advised in the",
				"use of scattering values obtained with transmittances between", 
				"0.3 and 0.1 and should be discarded below 0.1"), call. = FALSE)


}
