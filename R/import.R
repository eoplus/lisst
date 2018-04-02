
#' Read LISST data
#'
#' Read LISST data (processed or binary files) for the registered instruments of
#' supported models.
#'
#' @param fl    Path to processed or binary file (e.g., *.DAT in the 
#'              LISST-100(X)).
#' @param zscat Path to ASCII background data file for LISST-100(X) (e.g., 
#'              *.asc). Ignored for LISST-200X. Only required when reading 
#'              binary files.
#' @param pl    Path length in meters. If not provided the function will assume 
#'              standard path length for the instrument model (with a warning), 
#'              i.e., no path reduction module in use.
#' @param sn    Serial number of the instrument. Optionally, can be omitted when 
#'              reading a SOP processed file or when out is 'raw' or 'cor' for 
#'              binary files. In this case, instrument model must be supplied to 
#'              argument model.
#' @param yr    The year of first measurement in the data file. Only necessary 
#'              for LISST-100(X). If not provided, function will make a best 
#'              guess (with a warning) based on the curent date and the julian
#'              date of first measurement.
#' @param out   The output format. Valid options are 'vol', 'pnc', 'raw', 'cor' 
#'              and 'cal'. See details.
#' @param model A character vector of the instrument model (e.g. "200" for 
#'              LISST-200X). For the LISST-100(X), the detector type must be 
#'              included in the name (e.g., "100C" or "100CX"). Ignored if sn is 
#'              provided.
#'
#' @details The function will determine the file type based on its extenssion. 
#' Processed files created from LISST SOP for the LISST-100(X) have extension 
#' .asc and from LISST-200X, extension .csv. Binary files have extension .DAT 
#' and RBN for LISST-100(X) and LISST-200X, respectivelly. The extension is also
#' used to determine 
#'
#' The level of processing will depend on the output requested by the 
#' user. For processed LISST files as creayed by the LISST SOP, the volume 
#' concentration ('vol', ppm) or the particle number concentration ('pnc', 
#' 1/L/µm) can be returned. For binary files, the raw digital counts ('raw'), 
#' the corrected digital counts ('cor') or the calibrated values ('cal') can be
#' returned. Corrected digital counts are the digital counts of the ring 
#' detectors de-attenuated, background corrected and compensated for ring area 
#' deviations from ideal log increase behaivour. All other parameters in 'cor' 
#' are kept at the original digital counts. Finally, for 'cal', all parameters 
#' are converted to physical units using the calibration constants registered 
#' with \code{lisst_reg} for the instrument serial number. In addition, types 
#' 'cal','vol' and 'pnc' return an extra column, with dates/times converted to 
#' POSIXct format.
#'
#' To allow easier manipulation, conversion between types, and plotting, all 
#' associated information for the instrument and model are saved as attributes 
#' in the returned object.
#'
#' @return A lisst object with attribute type according to the processing level 
#' requested. 
#'
#' @references
#' MATLAB source code provided by Sequoia Scientific, Inc, and available at:
#' https://www.sequoiasci.com/product/lisst-100x/, 
#' https://www.sequoiasci.com/product/lisst-200x/
#'
#' @export

read_lisst <- function(fl, sn, pl, zscat, yr, out, model) {
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

	if(missing(sn)) {
		if(out == "cor" || out == "cal")
			stop(paste("Processing of binary to 'cor' or 'cal' requires instrument specific", 
				"information. sn must be provided."), call. = FALSE)
		else if(missing(model))
			stop(paste("For reading processed files or for 'raw' and 'cor' outputs",
				"from binary files, model must be supplied if sn is not"), call. = FALSE)

		linst <- list(X = FALSE)
		if(length(grep("X", model)) > 0)
			linst <- list(X = TRUE)
		model <- sub("X", "", model)
		lmodl <- switch(model, 
			"200"  = .getmodp(list(mod = "200", dty = "A")),
			"100B" = .getmodp(list(mod = "100", dty = "B")),
			"100C" = .getmodp(list(mod = "100", dty = "C"))
		)
		if(is.null(lmodl))
			stop("model must be one of '100B(X)', '100C(X)' or '200'", call. = FALSE)
		sn    <- NA
	} else {
		sn <- as.character(sn)
		linst <- .LISSTi[[sn]]
		if(is.null(linst))
			stop("Instrument not registered. See ?lisst_reg", call. = FALSE)
		lmodl <- .getmodp(linst)
		model <- NULL
	}
	if(missing(zscat)) {
		if(out == 'cor' || out == 'cal')
			stop("zscat must be provided for out 'cor' or 'cal'", call. = FALSE)
		else
			zscat <- NULL
	} else if(!file.exists(zscat)) {
		if(out == 'cor' || out == 'cal')
			stop(paste("File", zscat, "not found"), call. = FALSE)
		else {
			zscat <- NULL
			warning(paste("File", zscat, "not found; zscat data will not be",
				"added to lisst object"), call. = FALSE)
		}
	}

	if(missing(pl)) {
		warning(paste("pl not provided; assuming standard path length of instrument model"), 
			call. = FALSE)
		pl <- lmodl$pl
	} else if(units::set_units(pl, m) > lmodl$pl) {
		stop(paste0("Path length in LISST-", linst$mod, " cannot be larger than ", 
			lmodl$pl, " m"), call. = FALSE)
	}

	if(missing(yr)) yr <- NULL

	if(mode == "binary") {
		if(lmodl$mod == "100" && grep('.\\.DAT', fl, perl = TRUE) < 1)
			stop("LISST-100(X) binary files must have a .DAT extension", call. = FALSE)
		if(lmodl$mod == "200" && grep('.\\.RBN', fl, perl = TRUE) < 1)
			stop("LISST-200X binary files must have a .RBN extension", call. = FALSE)
		lo <- lisst_bin(fl = fl, sn = sn, pl = pl, zscat = zscat, linst = linst, lmodl = lmodl)
		if(out == 'cal') lo <- lgetcal(lo, yr)
		else if(out == 'cor') lo <- lgetcor(lo)
	} else {
		lo <- lisst_pro(fl = fl, sn = sn, pl = pl, zscat = zscat, linst = linst, lmodl = lmodl)
		if(out == 'pnc') lo <- lgetpnc(lo)
	}

	return(lo)
}

#' Read a LISST processed file
#'
#' Read a LISST processed file for the registered instruments of supported 
#' models. It is not intended to be used directly, but called from 
#' \code{read_lisst}.

lisst_pro <- function(fl, sn, pl, zscat, linst, lmodl, yr) {

	if(lmodl$mod == "100") {
                if(length(grep("_rs", fl)) > 0) ity <- "rs"
		else ity <- "ss"
		lo <- read.table(fl, header = FALSE)
	}
	if(lmodl$mod == "200") {
                if(length(grep("_rs", fl)) > 0) ity <- "rs"
		else ity <- "ss"		
		lo <- read.csv(fl, header = FALSE)
	}
	lo <- as.data.frame(lo)
	colnames(lo) <- lmodl$lvarn
	lo[, "Beam attenuation"] <- lo[, "Beam attenuation"] * as.numeric(lmodl$pl / pl)
	for(i in c(1:38, 41:42)) units(lo[, i]) <- lmodl$varun[i]

	zscatd <- rep(NA, lmodl$bnvar)
	if(!(missing(zscat) || is.null(zscat)))
		if(file.exists(zscat))
			zscatd <- as.numeric(read.table(zscat)[, 1])
	names(zscatd) <- lmodl$lvarn[1:lmodl$bnvar]

	linst$ity <- ity
	attr(lo, "sn")    <- sn
	attr(lo, "type")  <- "vol"
	attr(lo, "linst") <- linst
	attr(lo, "lmodl") <- lmodl
	attr(lo, "zscat") <- zscatd
	attr(lo, "class") <- c("lisst", "data.frame")
	lo$Time <- lgdate(lo, yr)
	return(lo)
}

#' Read a LISST binary file
#'
#' Read a LISST binary file for the registered instruments of supported models. 
#' It is not intended to be used directly, but called from \code{read_lisst}.

lisst_bin <- function(fl, sn, pl, zscat, linst, lmodl) {
	if(lmodl$mod == "100") {
		lo <- readBin(fl, "integer", n = file.info(fl)$size, size = 1, signed = FALSE, 
			endian = "little")
		lo <- lo[seq(1, length(lo), 2)] * 256 + lo[seq(2, length(lo), 2)]
		nrows <- floor(length(lo) / lmodl$bnvar)
		lo <- lo[1:(nrows * lmodl$bnvar)]
		lo <- matrix(lo, nrow = nrows, byrow = TRUE)
		if(linst$X) {
			lo[, 1:32] <- lo[, 1:32] / 10
		}

		zscatd <- rep(NA, lmodl$bnvar)
		if(!(missing(zscat) || is.null(zscat)))
			if(file.exists(zscat))
				zscatd <- as.numeric(read.table(zscat)[, 1])
		names(zscatd) <- lmodl$lvarn[1:lmodl$bnvar]

		lo <- as.data.frame(lo)
		colnames(lo) <- lmodl$lvarn[1:lmodl$bnvar]
		attr(lo, "sn")    <- sn
		attr(lo, "type")  <- "raw"
		attr(lo, "linst") <- linst
		attr(lo, "lmodl") <- lmodl
		attr(lo, "zscat") <- zscatd
		attr(lo, "class") <- c("lisst", "data.frame")

		return(lo)
	}

	if(ancd$mod == "200") {
		cat("To be added soon...")
	}
}


