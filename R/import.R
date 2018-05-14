
#' Read LISST data
#'
#' Read LISST processed or binary data files. 
#'
#' @param fl    Path to processed or binary file (e.g., *.DAT for the 
#'              LISST-100(X)).
#' @param zscat Background data. Might be provided as a raw lisst object or as a 
#'              path to file (e.g., *.asc for LISST-100(X)). This is required 
#'		for processing binary files, with the exception of out = 'raw'. 
#'		It is optional for the LISST-200X, if it is desired to process 
#'		the binary file with a different zscat than the internaly 
#'              stored.
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
#' @param out   The output format. Valid options are 'raw', 'cor', 'cal', 'vsf', 
#'              'vol' and 'pnc'. See details.
#' @param model A character vector of the instrument model (e.g. "200" for 
#'              LISST-200X). For the LISST-100(X), the detector type must be 
#'              included in the name (e.g., "100C" or "100XC"). Ignored if sn is 
#'              provided.
#' @param tz    Time zone to interpret time indexing. Defaults to UTC.
#' @param trant Logical. Should a optical transmittance threshold be applied? If
#'              TRUE, ring values for a given sample with transmittance < 0.3 
#'              are set to NA.
#'
#' @details
#' The function will determine the file type based on its extenssion. Processed 
#' files created from LISST SOP for the LISST-100(X) have extension .asc and 
#' for the LISST-200X, extension .csv. Binary files have extension .DAT and .RBN 
#' for the LISST-100(X) and LISST-200X, respectivelly.
#'
#' The parameter out determine the level of processing of the returned object. 
#' For binary files out can be 'raw' for the raw digital counts, 'cor' for the 
#' corrected digital counts, 'cal' for values in calibrated physical units or 
#' 'vsf' for the volume scattering function. The corrected digital counts are 
#' the raw counts de-attenuated for the particle \strong{and water} extinction, 
#' background subtracted and compensated for area deviations from nominal 
#' values. 'cal' applies the instrument specific calibration constants to 'cor' 
#' (for all variables). Aditionally, the transmittance and the particle beam 
#' attenuation are added in 'cal' type lisst objects. The particle volume 
#' scattering function is the calibrated values normalized by the detector solid
#' angle, length of the path generating the signal and energy enetering the 
#' path. If out is not provided, 'vsf' will be returned for a binary file input. 
#' For processed files, out can be 'vol' for the volume concentration (ppm) or 
#' 'pnc' for the number concentration (1/L/µm). If not provided, 'vol' will be 
#' returned for processed files. As of this version, is not possible to direcly 
#' retrieve the particle size distribution (PSD) from binary data (inversion 
#' model not implemented), so 'vol' and 'pnc' can only be selected for processed 
#' files. Functions \code{lget}, \code{lgetraw}, \code{lgetcor}, \code{lgetcal}, 
#' \code{lgetvsf}, \code{lgetvol} and \code{lgetpnc} allow to switch between 
#' types without need to read from disk.
#'
#' If the background data is saved as a binary file, it must be opened first 
#' with \code{read_lisst} and aggregated to one row lisst object with 
#' \link{lstat} before being passed to a new \code{read_lisst} call to open the 
#' data itself. 
#'
#' The time indexing, with date/time in POSIX format is added to all created 
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
#' @return 
#' The function returns an object of class lisst, an S3 Class inheriting from 
#' classes data.frame and units. Essentially the data is stored as units objects 
#' in a data.frame, wich stores the ancilary data as attributes. The number of 
#' columns in the data.frame is model dependent and is always one larger than 
#' the standard LISST data since date/time in POSIXct format is added to the 
#' objects.
#'
#' The attributes are not expected to be manipulated directly, so are only 
#' briefly described:
#' \itemize{
#'   \item type  - The type of data in the lisst object;
#'   \item linst - Instrument specific data;
#'   \item lmodl - Model specific data;
#'   \item lproc - Processing applied to data (including inversion model);
#'   \item zscat - Background scattering values.
#' }
#'
#' @references
#' MATLAB source code provided by Sequoia Scientific, Inc, and available at:
#' https://www.sequoiasci.com/product/lisst-100x/, 
#' https://www.sequoiasci.com/product/lisst-200x/
#'
#' @seealso \code{\link{lget}}, \code{\link{lgetraw}}, \code{\link{lgetcor}}, 
#' \code{\link{lgetcal}}, \code{\link{lgetvsf}}, \code{\link{lgetvol}} and 
#' \code{\link{lgetpnc}}.
#'
#' @examples
#' flp <- system.file("extdata", "DN_27_rs.asc", package = "lisst")
#' flb <- system.file("extdata", "DN_27.DAT", package = "lisst")
#'
#' # For a unregistered LISST instrument:
#' model <- "100CX" 
#' lop <- read_lisst(flp, model = model)
#' lob <- read_lisst(flb, out = 'raw', model = model)
#' 
#' # If other levels of processing, including VSF, are required, first
#' # register a LISST instrument:
#' path  <- system.file("extdata", package = "lisst")
#' model <- 100
#' lisst_reg(model, path)
#'
#' sn    <- 1298
#' pl    <- 0.05
#' yr    <- 2018
#' out   <- 'vol'
#' zscat <- system.file("extdata", "bg_20180326.asc", package = "lisst")
#'
#' # For a processed file:
#' lop <- read_lisst(flp, sn, pl, zscat, yr, out)
#' lop <- read_lisst(flp, sn, pl, zscat, yr)
#' lop <- read_lisst(flp, sn, pl, zscat)
#' lop <- read_lisst(flp, sn, pl)
#' lop <- read_lisst(flp, sn) # minimum information
#'
#' # For a binary file:
#' lob <- read_lisst(flb, sn, pl, zscat, yr, out = 'raw')
#' lob <- read_lisst(flb, sn, pl, zscat, yr, out = 'cor')
#' lob <- read_lisst(flb, sn, pl, zscat, yr, out = 'cal')
#' lob <- read_lisst(flb, sn, pl, zscat, yr)
#' lob <- read_lisst(flb, sn, pl, zscat)
#' lob <- read_lisst(flb, sn, zscat = zscat) # minimum information for full capability
#'
#' @export

read_lisst <- function(fl, sn, pl, zscat, yr, out, model, tz = 'UTC', trant = TRUE) {
	if(!file.exists(fl))
		stop(paste("File", fl, "not found"), call. = FALSE)

	mode <- "processed"
	if(length(grep('.(\\.asc|\\.csv)', fl, perl = TRUE)) < 1)		# Find if file is processed or binary by extension
		mode <- "binary"

	if(missing(out) && mode == "processed") out <- "vol"			# Default output for processed files
	if(missing(out) && mode == "binary")    out <- "vsf"			# Default output for binary files
	if(mode == "processed" && !(out == "vol" || out == "pnc"))
		stop("out for LISST SOP processed files must be 'vol' or 'pnc'", 
			call. = FALSE)
	if(mode == "binary" && (out == "vol" || out == "pnc"))
		stop("out for binary LISST files must be 'raw', 'cor','cal' or", 
			" vsf", call. = FALSE)

	# Retrieve information on model, either from sn or model parameters:
	if(missing(sn)) {
		if(out == "cor" || out == "cal" || out == 'vsf')
			stop("Processing of binary to 'cor', 'cal' or 'vsf' ",
				"requires instrument specific information. sn ",
				"of a registered instrument must be provided.", 
				call. = FALSE)
		else if(missing(model))
			stop("For reading processed files or for 'raw' outputs", 
				" from binary files, model must be supplied if", 
				" sn is not", call. = FALSE)

		linst <- list(X = FALSE)
		if(length(grep("X", model)) > 0) linst$X <- TRUE
		model <- sub("X", "", model)
		linst$mod <- model
		lmodl <- switch(model,
			"200"  = .getmodp(list(mod = "200", dty = "A")),
			"100B" = .getmodp(list(mod = "100", dty = "B")),
			"100C" = .getmodp(list(mod = "100", dty = "C"))
		)
		if(is.null(lmodl))
			stop("model must be one of '100(X)B', '100(X)C' or ", 
				"'200X'", call. = FALSE)
		sn    <- NA
	} else {
		sn <- as.character(sn)
		linst <- .LISSTi[[sn]]
		if(is.null(linst))
			stop("Instrument not registered. See ?lisst_reg", 
				call. = FALSE)
		lmodl <- .getmodp(linst)
		model <- NULL
	}

	# Check that zscat is necessary and if exists:
	if(missing(zscat)) {
		if((out == 'cor' || out == 'cal' || out == 'vsf') && linst$mod == '100')
			stop("zscat must be provided for out 'cor', 'cal' or ", 
				"'vsf'", call. = FALSE)
		else
			zscat <- NULL
	} else if(is.character(zscat)) {
		if(!file.exists(zscat)) {
			if(out == 'cor' || out == 'cal' || out == 'vsf')
				stop(paste("File", zscat, "not found"), 
					call. = FALSE)
			else {
				zscat <- NULL
				warning(paste("File", zscat, "not found; zscat", 
					" data will not be added to lisst ", 
					"object"), call. = FALSE)
			}
		}
	}

	# lmodl$pl is not set directly here from pl parameter because .lisst_pro
	# will need to compensate the beam attenuation (and the volume 
	# concentration) with a factor based on the standard pl of the model:
	if(missing(pl)) {
		warning(paste("pl not provided - assuming standard path length"), 
			call. = FALSE)
		pl <- lmodl$pl
	} else if(pl > drop_units(lmodl$pl)) {
		stop(paste0("Path length in LISST-", linst$mod, " cannot be ", 
			"larger than ", drop_units(lmodl$pl), " m"), 
			call. = FALSE)
	}
	pl <- set_units(pl, 'm')

	# Manage year for 100(X) models:
	guess <- FALSE
	if(missing(yr)) {
		if(lmodl$mod == "100")
			warning("yr not provided - using best guess", 
				call. = FALSE)
		yr  <- format(file.info(fl)$mtime, "%Y")
		guess <- TRUE
	}

	# Call specific import functions:
	if(mode == "binary") {
		if(lmodl$mod == "100" && grep('.\\.DAT', fl, perl = TRUE) < 1)
			stop("LISST-100(X) binary files must have a .DAT ", 
				"extension", call. = FALSE)
		if(lmodl$mod == "200" && grep('.\\.RBN', fl, perl = TRUE) < 1)
			stop("LISST-200X binary files must have a .RBN ", 
				"extension", call. = FALSE)
		lo <- .lisst_bin(fl = fl, sn = sn, pl = pl, zscat = zscat, 
			linst = linst, lmodl = lmodl)
	} else {
		if(lmodl$mod == "100" && grep('.\\.asc', fl, perl = TRUE) < 1)
			stop("LISST-100(X) processed files must have a .asc ", 
				"extension", call. = FALSE)
		if(lmodl$mod == "200" && grep('.\\.csv', fl, perl = TRUE) < 1)
			stop("LISST-200X processed files must have a .csv ", 
				"extension", call. = FALSE)
		lo <- .lisst_pro(fl = fl, sn = sn, pl = pl, zscat = zscat, 
			linst = linst, lmodl = lmodl)
	}

	# Get time reference:
	ti <- .lgdate(lo, yr, tz, guess)
	tr <- as.numeric(difftime(max(ti), min(ti), units = "sec"))
	ti <- rep(ti[1], nrow(lo)) + 						# Add precision to avoid duplicate rownames
		cumsum(c(0, rep(tr / (nrow(lo) - 1), nrow(lo) - 1)))

	rownames(lo) <- format(ti, "%Y-%m-%d %H:%M:%OS1 %Z")
	lo <- lget(lo, out)
	lo <- lo[, -which(colnames(lo) %in% c("Time1", "Time2", "Year", 	# Remove the original time columns and not used column in LISST-200X
		"Month", "Day", "Hour", "Minute", "Second", "NU"))]

	# Apply simple quality check:
	if(mode == 'binary' & is.null(zscat)) trant <- FALSE
	if(trant) {
		if(out == 'raw' | out == 'cor')
			ot <- lget(lo, 'cal')[, "OptTrans"]
		else
			ot <- lo[, "OptTrans"]
		id  <- which(drop_units(ot) < 0.3)
		for(i in id) lo[i, 1:attr(lo, 'lmodl')$nring] <- NA
	}

	return(lo)
}

# Read a LISST processed file
#
# Read a LISST processed file for the registered instruments of supported 
# models. It is not intended to be used directly, but called from 
# \code{read_lisst}.

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
	mfact <- drop_units(lmodl$pl / pl)					# PRM correction factor
	for(i in 1:lmodl$nring) lo[, i] <- lo[, i] * mfact			# PRM correction for volume concentration
	lo[, "BeamAtt"] <- lo[, "BeamAtt"] * mfact				# PRM correction for beam attenuation
	lmodl$pl <- pl								# Store the true path length in object metadata

	# Process background data:
	zscatd <- rep(as.numeric(NA), lmodl$bnvar)
	if(!(missing(zscat) || is.null(zscat)))
		if(is.lisst(zscat)) {
			if(nrow(zscat) > 1)
				stop("A lisst object used as zscat must have a", 
					" single row. Use lstat for ", 
					"aggregation.", call. = FALSE)
			zscatd <- drop_lisst(lget(zscat, 'raw'))
		} else if(file.exists(zscat)) {
			zscatd <- as.numeric(read.table(zscat)[, 1])
		}
	zscatd <- as.data.frame(matrix(zscatd, nrow = 1))
	names(zscatd) <- lmodl$lvarn[1:lmodl$bnvar]

	# Add units, avoiding data columns related to time:
	id <- setdiff(lmodl$lvarn, c("Time1", "Time2", "Year", "Month", "Day", 
		"Hour", "Minute", "Second"))
	id <- which(lmodl$lvarn %in% id)
	for(i in id) {
		units(lo[, i]) <- lmodl$varun[i]
		if(i < lmodl$bnvar) units(zscatd[, i]) <- 1
	}

	lo <- structure(lo, type = 'vol', lproc = list(ity = ity), 
		linst = linst, lmodl = lmodl, zscat = zscatd, 
		class = c("lisst", "data.frame"))
	return(lo)
}

# Read a LISST binary file
#
# Read a LISST binary file for the registered instruments of supported models. 
# It is not intended to be used directly, but called from \code{read_lisst}.
#
# Those are based on the Matlab functions provided by Sequoia.

.lisst_bin <- function(fl, sn, pl, zscat, linst, lmodl) {
	lmodl$pl <- pl
	if(lmodl$mod == "100") {
		lo <- readBin(fl, "integer", n = file.info(fl)$size, size = 1, 
			signed = FALSE, endian = "little")
		lo <- lo[seq(1, length(lo), 2)] * 256 + lo[seq(2,length(lo), 2)]
		nrows <- floor(length(lo) / lmodl$bnvar)
		lo <- lo[1:(nrows * lmodl$bnvar)]
		lo <- matrix(lo, nrow = nrows, byrow = TRUE)
		if(linst$X) {
			lo[, 1:32] <- lo[, 1:32] / 10
		}

		zscatd <- rep(as.numeric(NA), lmodl$bnvar)
		if(!(missing(zscat) || is.null(zscat)))
			if(is.lisst(zscat)) {
				if(nrow(zscat) > 1)
					stop("A lisst object used as zscat ", 
					"must have a single row. Use lstat for", 
					" aggregation.", call. = FALSE)
				zscatd <- drop_lisst(lget(zscat, 'raw'))
			} else if(file.exists(zscat)) {
				zscatd <- as.numeric(read.table(zscat)[, 1])
			}
		zscatd <- as.data.frame(matrix(zscatd, nrow = 1))
		names(zscatd) <- lmodl$lvarn[1:lmodl$bnvar]

		lo <- as.data.frame(lo)
		colnames(lo) <- lmodl$lvarn[1:lmodl$bnvar]
		id <- setdiff(lmodl$lvarn[1:lmodl$bnvar], c("Time1", "Time2", 
			"Year", "Month", "Day", "Hour", "Minute", "Second"))
		id <- which(lmodl$lvarn %in% id)
		for(i in id) {
			units(lo[, i]) <- 1
			units(zscatd[, i]) <- 1
		}

		lo <- structure(lo, type = 'raw', lproc = list(ity = NA), 
			linst = linst, lmodl = lmodl, zscat = zscatd, 
			class = c("lisst", "data.frame"))
		return(lo)
	}

	if(lmodl$mod == "200") {
		RecordSize <- 120						# 120 bytes per record
		num16 = (RecordSize / 2) - 1
		num32 = floor((RecordSize - 2) / 4)

		fid  <- file(fl, "rb")
		seek(fid, (8 * RecordSize) + 2, 'start')
		zsc <- readBin(fid, "integer", n = num16, size = 2, signed = FALSE, endian = "big")
		seek(fid, (9 * RecordSize) + 2, 'start')
		lo <- readBin(fid, "integer", n = file.info(fl)$size - seek(fid), size = 2, signed = FALSE, endian = "big")
		lo <- matrix(c(lo, 56026), ncol = 60, byrow = T)[, -60]
		close(fid)

		lo[lo > 40950] <- lo[lo > 40950] - 65536			# Check this, why? and why 65536 instead of 65535?
		zsc[zsc > 40950] <- zsc[zsc > 40950] - 65536
		lo[, 1:lmodl$nring] <- lo[, 1:lmodl$nring] / 10 
		zsc[1:lmodl$nring] <- zsc[1:lmodl$nring] / 10

		zscatd <- rep(as.numeric(NA), lmodl$bnvar)
		if(!(missing(zscat) || is.null(zscat))) {
			if(is.lisst(zscat)) {
				if(nrow(zscat) > 1)
					stop("A lisst object used as zscat ", 
					"must have a single row. Use lstat for", 
					" aggregation.", call. = FALSE)
				zscatd <- drop_lisst(lget(zscat, 'raw'))
			} else if(file.exists(zscat)) {
				zscatd <- as.numeric(read.table(zscat)[, 1])
			}
		} else {
			zscatd <- zsc
		}
		zscatd <- as.data.frame(matrix(zscatd, nrow = 1))
		names(zscatd) <- lmodl$lvarn[1:lmodl$bnvar]

		lo <- as.data.frame(lo)
		colnames(lo) <- lmodl$lvarn[1:lmodl$bnvar]
		id <- setdiff(lmodl$lvarn[1:lmodl$bnvar], c("Time1", "Time2", 
			"Year", "Month", "Day", "Hour", "Minute", "Second"))
		id <- which(lmodl$lvarn %in% id)
		for(i in id) {
			units(lo[, i]) <- 1
			units(zscatd[, i]) <- 1
		}

		lo <- structure(lo, type = 'raw', lproc = list(ity = NA), 
			linst = linst, lmodl = lmodl, zscat = zscatd, 
			class = c("lisst", "data.frame"))
		return(lo)
	}
}

# Get LISST measurement dates
# 
# The function retrieves the date/time of measurements for the records in the 
# lisst object as POSIXct objects. The time zone of the dates will depend on 
# user tz options. It is not intended to be used directly, but called from
# read_lisst function.
#
# param lo    A lisst object.
# param yr    The year of first measurement in the lisst object. Ignored for 
#             LISST-200X.
# param guess Logical. Is the yr passed a guess from read_lisst?
#
# details
# If yr was not provided by the user, read_lisst will call lgdate with guess = 
# TRUE. This informs the function that the year passed is an estimate of the 
# last year of measurement, not the first as the yr passed by the user. It was 
# kept this way since LISST users might be used to that way of informing date, 
# but the guess date has to be about the last measurement. See the details of 
# the read_lisst function on how the guess is made.

.lgdate <- function(lo, yr, tz, guess = FALSE) {
	lmodl <- attr(lo, "lmodl")
	lo    <- drop_lisst(lo) 
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
        	dates  <- as.POSIXct(paste(yr, julian, hour, min, sec, 
			sep = "-"), format = "%Y-%j-%H-%M-%S", tz = tz)
	} else if(lmodl$mod == "200") {
		dates  <- as.POSIXct(paste(lo[, 43], lo[, 44], lo[, 45], 
				lo[, 46], lo[, 47], lo[, 48], sep = "-"), 
				format = "%Y-%m-%d-%H-%M-%S", tz = tz)
	}
	return(dates)
}


