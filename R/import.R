#' @import units
NULL

#' Read LISST data
#'
#' Read LISST processed or binary data files. 
#'
#' @param fl    Path to processed or binary file (e.g., *.DAT for the 
#'              LISST-100(X)).
#' @param zscat Path to ASCII background data file (e.g., *.asc). This is 
#'              required for processing binary files, with the exception of out 
#'		= 'raw'. It is optional for the LISST-200X, if it is desired to  
#'              process the binary file with a different zscat than the 
#'              internaly stored.
#' @param pl    Path length in meters. If not provided the function will assume 
#'              standard path length for the instrument model (with a warning), 
#'              i.e., no path reduction module in use.
#' @param sn    Serial number of the instrument. Optionally, can be omitted when 
#'              reading a LISST-SOP processed file or when out = 'raw' (binary 
#'              files). In this case, instrument model must be supplied to 
#'              argument model.
#' @param yr    The year of first measurement in the data file. Only necessary 
#'              for LISST-100(X). If not provided, function will make a best 
#'              guess (with a warning). See details.
#' @param out   The output format. Valid options are 'raw', 'cor', 'cal', 'vol' 
#'              and 'pnc'. See details.
#' @param model A character vector of the instrument model (e.g. "200" for 
#'              LISST-200X). For the LISST-100(X), the detector type must be 
#'              included in the name (e.g., "100C" or "100XC"). Ignored if sn is 
#'              provided.
#'
#' @details
#' The function will determine the file type based on its extenssion. Processed 
#' files created from LISST SOP for the LISST-100(X) have extension .asc and 
#' for the LISST-200X, extension .csv. Binary files have extension .DAT and .RBN 
#' for the LISST-100(X) and LISST-200X, respectivelly.
#'
#' The parameter out determine the level of processing of the returned object. 
#' For binary files out can be 'raw' for the raw digital counts, 'cor' for the 
#' corrected digital counts or 'cal' for values in calibrated physical units. 
#' The corrected digital counts are the raw counts de-attenuated for the 
#' particle extinction, background subtracted and compensated for area 
#' deviations from nominal values. 'cal' applies the instrument specific 
#' calibration constants to 'cor' (for all variables). If out is not provided, 
#' 'cal' will be returned for a binary file input. For processed files, out can 
#' be 'vol' for the volume concentration (ppm) or 'pnc' for the number 
#' concentration (1/L/µm). If not provided, 'vol' will be returned for processed
#' files. As of this version, is not possible to direcly retrieve the particle 
#' size distribution (PSD) from binary data (inversion model not implemented), 
#' so 'vol' and 'pnc' can only be selected for processed files. Functions 
#' \code{lgetraw}, \code{lgetcor}, \code{lgetcal}, \code{lgetvol} and 
#' \code{lgetpnc} allow to switch between types without need to read from disk. 
#'
#' A column "Time", with date/time in POSIXct format is added to all created 
#' objects. If yr is missing when reading a LISST-100(X) file, the function will 
#' 'guess' its value, by acessing the file system modification date information. 
#' The modification date is used to be consistent across platforms, since 
#' UNIX-type systems do not register creation date. Still, those are expceted to 
#' be equivalent since the files are not expected to be modified since their 
#' creation by the LISST instrument (binary) or the LISST-SOP (processed). In 
#' the case of binary file, the year will be then precise for the year of 
#' \emph{last} measurement and the function will handle it appropriatly. For a 
#' processed file the logic will break down in cases that the file is processed 
#' in a different year than the final measurement. If a guess was not correct 
#' and the user have the approriate information, is possible to directly alter 
#' the Time column in the lisst object using standard R tools to handle POSIX 
#' objects as an alternative to call \code{read_lisst} again.
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

	# It is not setting directly the pl component of lmodl because .lisst_pro
	# will need to compensate the beam attenuation with a factor based on the 
	# standard pl of the model. 
	if(missing(pl)) {
		warning(paste("pl not provided - assuming standard path length"), 
			call. = FALSE)
		pl <- lmodl$pl
	} else if(units::set_units(pl, m) > lmodl$pl) {
		stop(paste0("Path length in LISST-", linst$mod, " cannot be larger than ", 
			lmodl$pl, " m"), call. = FALSE)
	}

	guess <- FALSE
	if(missing(yr)) {
		if(lmodl$mod == "100")
			warning("yr not provided - using best guess", call. = FALSE)
		yr  <- format(file.info(fl)$mtime, "%Y")
		guess <- TRUE
	}


	if(mode == "binary") {
		if(lmodl$mod == "100" && grep('.\\.DAT', fl, perl = TRUE) < 1)
			stop("LISST-100(X) binary files must have a .DAT extension", call. = FALSE)
		if(lmodl$mod == "200" && grep('.\\.RBN', fl, perl = TRUE) < 1)
			stop("LISST-200X binary files must have a .RBN extension", call. = FALSE)
		lo <- .lisst_bin(fl = fl, sn = sn, pl = pl, zscat = zscat, linst = linst, lmodl = lmodl)
		if(out == 'cal') lo <- lgetcal(lo)
		else if(out == 'cor') lo <- lgetcor(lo)
	} else {
		lo <- .lisst_pro(fl = fl, sn = sn, pl = pl, zscat = zscat, linst = linst, lmodl = lmodl)
		if(out == 'pnc') lo <- lgetpnc(lo)
	}

	lo$Time <- .lgdate(lo, yr, guess)
	return(lo)
}

#' Read a LISST processed file
#'
#' Read a LISST processed file for the registered instruments of supported 
#' models. It is not intended to be used directly, but called from 
#' \code{read_lisst}.

.lisst_pro <- function(fl, sn, pl, zscat, linst, lmodl) {

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
	lmodl$pl <- pl

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
	return(lo)
}

#' Read a LISST binary file
#'
#' Read a LISST binary file for the registered instruments of supported models. 
#' It is not intended to be used directly, but called from \code{read_lisst}.

.lisst_bin <- function(fl, sn, pl, zscat, linst, lmodl) {
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

#' Get LISST measurement dates
#' 
#' The function retrieves the date/time of measurements for the records in the 
#' lisst object as POSIXct objects. The time zone of the dates will depend on 
#' user tz options. It is not intended to be used directly, but called from
#' read_lisst function.
#'
#' @param lo    A lisst object.
#' @param yr    The year of first measurement in the lisst object. Ignored for 
#'              LISST-200X.
#' @param guess Logical. Is the yr passed a guess from read_lisst?
#'
#' @details
#' If yr was not provided by the user, read_lisst will call lgdate with guess = 
#' TRUE. This informs the function that the year passed is an estimate of the 
#' last year of measurement, not the first as the yr passed by the user. It was kept
#' this way since LISST users might be used to that way of informing date, but the
#' guess date has to be about the last measurement. See the details of the 
#' read_lisst function on how the guess is made.

.lgdate <- function(lo, yr, guess = FALSE) {
	lmodl <- attr(lo, "lmodl")
	if(lmodl$mod == "100") {
		julian <- floor(lo[, 39] / 100)
		id <- which(diff(julian) <= -364)
        	yr <- rep(yr, nrow(lo))
		if(length(id) > 0) {
			if(guess) yr <- yr - length(id)
			for(i in 1:length(id)) {
				yr[id[i]:length(yr)] <- yr[id[i]:length(yr)] + 1
			}
		}
		hour   <- round(((lo[, 39] / 100) - julian) * 100)
		min    <- floor(lo[, 40] / 100)
        	sec    <- round(((lo[, 40] / 100) - min) * 100)
        	dates  <- as.POSIXct(paste(yr, julian, hour, min, sec, sep = "-"), 
			format = "%Y-%j-%H-%M-%S")
	} else if(lmodl$mod == "200") {
		dates  <- as.POSIXct(paste(lo[, 43], lo[, 44], lo[, 45], lo[, 46], lo[, 47], 
				lo[, 48], sep = "-"), format = "%Y-%m-%d-%H-%M-%S")
	}
	return(dates)
}

