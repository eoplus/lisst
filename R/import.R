
#' Read LISST data
#'
#' Read LISST data (processed or binary files) for the registered instruments of
#' supported models.
#'
#' @param fl    Path to binary file (e.g., *.DAT in the LISST-100(X)).
#' @param zscat Path to ASCII background data file for LISST-100(X) (e.g., 
#'              *.asc). Ignored for LISST-200X.
#' @param pl    Path length in meters. If not provided the system will assumed 
#'              standard path length for the instrument model, i.e., no path 
#'              reduction module in use.
#' @param sn    Serial number of the instrument.
#' @param yr    The year of first measurement in the data file. Only necessary 
#'              for LISST-100(X).
#' @param out   The output format. Valid options are 'vol', 'psd', 'raw', 'cor' 
#'              and 'cal'. See details.
#'
#' @details The function will determine the file type based on its extenssion. 
#' Processed files outputed from LISST SOP for the LISST-100(X) have extension 
#' .asc and from LISST-200X, extension .csv. Binary files have extension .DAT 
#' and RBN for LISST-100(X) and LISST-200X, respectivelly.
#'
#' The level of processing will depend on the output requested by the 
#' user. For processed LISST files as outputed by the LISST SOP the volume 
#' concentration ('vol', ppm) or the particle number concentration ('pnc', 
#' 1/L/µm) can be returned. For binary files, the raw digital counts ('raw'), 
#' the corrected digital counts ('cor') or the calibrated values ('cal') can be
#' returned. Corrected digital counts are the digital counts of the ring 
#' detectors de-attenuated, background corrected and compensated for ring area 
#' deviations from ideal log increase behaivour. All other parameters are in 
#' 'cor' are kept at digital counts. Finally, for 'cal', all parameters are 
#' converted to physical units using the calibration constants registered with 
#' \code{lisst_reg} for the instrument serial number. 
#'
#' @return A lisst object with attribute type according to the processing level 
#' requested. 
#'
#' @export

read_lisst <- function(fl, zscat, sn, pl, yr, out) {
	sn <- as.character(sn)
	linst <- .LISSTi[[sn]]
	if(is.null(linst))
		stop("Instrument not registered. See ?lisst_reg", call. = FALSE)
	lmodl <- .getmodp(linst)

	if(!file.exists(fl))
		stop(paste("File", fl, "not found"), call. = FALSE)
	mode <- "processed"
	if(length(grep('.(\\.asc|\\.csv)', fl, perl = TRUE)) < 1) mode <- "binary"

	if(missing(out) && mode == "processed") out <- "vol"
	if(missing(out) && mode == "binary")    out <- "cal"
	if(mode == "processed" && !(out == "vol" || out == "pnc"))
		stop("out for processed SOP files must be 'vol' or 'pnc'", call. = FALSE)
	if(mode == "binary" && (out == "vol" || out == "pnc"))
		stop("out for binary LISST files must be 'raw', 'cor' or 'cal'", call. = FALSE)

	if(!file.exists(zscat) && mode == "binary")
		stop(paste("File", zscat, "not found"), call. = FALSE)
	if(missing(pl)) {
		warning(paste("pl not provided; assuming standard path length of instrument model"), 
			call. = FALSE)
		pl <- lmodl$pl
	} else {
		if(units::set_units(pl, m) > lmodl$pl)
			stop(paste0("Path length in LISST-", linst$mod, " cannot be larger than ", 
				lmodl$pl, " m"), call. = FALSE)
	}
	if(missing(yr)) yr <- NULL

	if(mode == "binary") {
		if(!file.exists(zscat))
			stop(paste("File", zscat, "not found"), call. = FALSE)
		if(linst$mod == "100" && grep('.\\.DAT', fl, perl = TRUE) < 1)
			stop("LISST-100(X) binary files must have a .DAT extension", call. = FALSE)
		if(linst$mod == "200" && grep('.\\.RBN', fl, perl = TRUE) < 1)
			stop("LISST-200X binary files must have a .RBN extension", call. = FALSE)
		lo <- lisst_bin(fl = fl, zscat = zscat, sn = sn, pl = pl, out = out)
	} else {
		if(!file.exists(zscat))
			zscat <- NULL
		lo <- lisst_pro(fl = fl, sn = sn, pl = pl, zscat = zscat)
	}

	typ <- attr(lo, "type")
	if(typ == "cal" || typ == "vol" || typ == "pnc") {
		lo$Time <- lgdate(lo, yr)
		if(linst$mod == "100") lo <- lo[, c(1:38, 41:43)]
		if(linst$mod == "200") lo <- lo[, c(1:42, 49:61)]
	}
	return(lo)
}

#' Read a LISST processed file
#'
#' Read a LISST processed file for the registered instruments of supported 
#' models. It is not intended to be used directly, but called from 
#' \code{read_lisst}.

lisst_pro <- function(fl, sn, pl, zscat) {
	linst <- .LISSTi[[sn]]
	lmodl <- .getmodp(linst)

	if(linst$mod == "100") lo <- read.table(fl, header = FALSE)
	if(linst$mod == "200") lo <- read.csv(fl, header = FALSE)

	lo <- as.data.frame(lo)
	colnames(lo) <- lmodl$lvarn
	lo[, "Beam attenuation"] <- lo[, "Beam attenuation"] * as.numeric(lmodl$pl / pl)
	for(i in c(1:38, 41:42)) units(lo[, i]) <- lmodl$varun[i]

	zscatd <- rep(NA, lmodl$bnvar)
	if(!(missing(zscat) || is.null(zscat)))
		if(file.exists(zscat))
			zscatd <- as.numeric(read.table(zscat)[, 1])

	attr(lo, "sn")    <- sn
	attr(lo, "type")  <- "vol"
	attr(lo, "linst") <- linst
	attr(lo, "lmodl") <- lmodl
	attr(lo, "zscat") <- zscatd
	attr(lo, "class") <- c("lisst", "data.frame")
	lo
}

#' Read a LISST binary file
#'
#' Read a LISST binary file for the registered instruments of supported models. 
#' It is not intended to be used directly, but called from \code{read_lisst}.

lisst_bin <- function(fl, zscat, sn, pl, out) {
	linst <- .LISSTi[[sn]]
	lmodl <- .getmodp(linst)

	if(linst$mod == "100") {
		lo <- readBin(fl, "integer", n = file.info(fl)$size, size = 1, signed = FALSE, 
			endian = "little")
		lo <- lo[seq(1, length(lo), 2)] * 256 + lo[seq(2, length(lo), 2)]
		nrows <- floor(length(lo) / lmodl$bnvar)
		lo <- lo[1:(nrows * lmodl$bnvar)]
		lo <- matrix(lo, nrow = nrows, byrow = TRUE)
		if(linst$X) {
			lo[, 1:32] <- lo[, 1:32] / 10
		}

		zscat <- as.numeric(read.table(zscat)[, 1])

		lo <- as.data.frame(lo)
		colnames(lo) <- lmodl$lvarn[1:lmodl$bnvar]
		attr(lo, "sn")    <- sn
		attr(lo, "type")  <- "raw"
		attr(lo, "linst") <- linst
		attr(lo, "lmodl") <- lmodl
		attr(lo, "zscat") <- zscat
		attr(lo, "class") <- c("lisst", "data.frame")

		if(out == "raw")
			return(lo)

		tau <- lo[, 33] * zscat[36] / zscat[33] / lo[, 36]
		lo[, 1:32] <- lo[, 1:32] / tau - (rep(zscat[1:32], each = nrow(lo)) * 
			lo[, 36] / zscat[36])
		lo[, 1:32] <- rep(linst$ringcf, each = nrow(lo)) * lo[, 1:32]
		lo[, 1:32][lo[, 1:32] < 0] <- 0

		attr(lo, "type")  <- "cor"

		if(out == "cor")
			return(lo)

		lo[, 33] <- lo[, 33] * linst$lpowcc[1] + linst$lpowcc[2]
		lo[, 34] <- lo[, 34] * linst$battcc[1] + linst$battcc[2]
		lo[, 35] <- lo[, 35] * linst$extrcc[1] + linst$extrcc[2]
		lo[, 36] <- lo[, 36] * linst$lrefcc[1] + linst$lrefcc[2]
		lo[, 37] <- lo[, 37] * linst$dpthcc[1] + linst$dpthcc[2]
		lo[, 38] <- lo[, 38] * linst$tempcc[1] + linst$tempcc[2]
		lo[lmodl$lvarn[41]] <- tau
		lo[lmodl$lvarn[42]] <- -log(tau) / lmodl$pl
		lo[, 1:32] <- as.data.frame(set_units(as.matrix(lo[, 1:32]) * linst$ringcc, mW))
		attr(lo, "type")  <- "cal"
		return(lo)
	}
	if(ancd$mod == "200") {
		cat("To be added soon...")
	}
}


