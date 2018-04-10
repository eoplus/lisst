
#' Register a new LISST instrument
#'
#' Register a new LISST instrument by importing calibration data and storing it 
#' for convenient processing of LISST data.
#'
#' @param model The model of the LISST instrument (e.g., "100" for the LISST-100
#'              and LISST-100X).
#' @param path  Path to directory containing instrument instalation files. See 
#'              details.
#'
#' @details For the LISST-100(X) the function expects that files follow the 
#' naming convention for LISST installation files. Specifically, it will search 
#' for a factory background file (factory_zsc_*.asc), a ring area file 
#' (ringarea_*.asc), an instrument data file (InstrumentData.txt) and a 
#' LISST-SOP ini file (lisst.ini). All files should be in the same directory 
#' provided in the path argument, and it should contain files for a single 
#' instrument.
#'
#' @examples:
#' path  <- system.file("extdata", package = "lisst")
#' model <- 100
#' lisst_reg(model, path)
#'
#' @export

lisst_reg <- function(model, path) {

	if(nargs() < 2)
		stop("lisst_reg has no defaults; all arguments must be specifyed", call. = FALSE)

	model <- sub("x", "", as.character(model))
	
	if(length(.LISSTm[[model]]) < 1)
		stop(paste("Model", model, "not supported. Supported models:", 
			paste(names(.LISSTm), collapse = ", ")), call. = FALSE)
	if(!dir.exists(path))
		stop(paste("Directory", path, "not found"), call. = FALSE)

	if(model == "100") {
		fls <- list.files(path, pattern = "", full.names = TRUE)
		zsc <- fls[grep("factory_zsc_", fls)]
		rig <- fls[grep("ringarea_", fls)]
		dat <- fls[grep("InstrumentData.txt", fls)]
		ini <- fls[grep("lisst.ini", fls)]

		if(length(zsc) == 0)
			stop("Factory background file not found", call. = FALSE)
		if(length(rig) == 0)
			stop("Ring area file not found", call. = FALSE)
		if(length(dat) == 0)
			stop("Instrument data file not found", call. = FALSE)
		if(length(ini) == 0)
			stop("LISST ini file not found", call. = FALSE)

		zsc <- read.table(zsc, header = FALSE)[, 1]
		rig <- as.numeric(read.table(rig, header = FALSE))
		dat <- read.csv(dat, header = FALSE)
		ini <- readLines(ini)

		hk  <- matrix(NA, ncol = 2, nrow = 10)
		rownames(hk) <- paste0("HK", 0:9)
		for(i in 0:9) {
			id <- grep(paste0("HK", i), ini)
			if(length(id) == 0) {
				hk <- hk[1:(i), ]
				break
			}
			rownames(hk)[i+1] <- strsplit(ini[id[1]], "=")[[1]][2] 
			hk[i+1, ] <- as.numeric(sapply(strsplit(ini[id[3:4]], "="), '[[', 2))
		}
		hkn <- rownames(hk)

		.LISSTi[[as.character(dat[1])]] <- list(
			mod = model,
			sn  = dat[, 1],
			dty = as.character(dat[, 2]),
			X   = dat[, 5] == "X",
			fzscat = zsc,
			ringcf = quantities::set_quantities(rig, 1, 0),
			ringcc = quantities::set_quantities(5 / 4096 / dat[, 3], W, 0),
			volcc  = dat[, 4],
			lpowcc = quantities::set_quantities(hk[grep("Laser Power", hkn), ], mW, 0),
			battcc = quantities::set_quantities(hk[grep("Battery", hkn), ], V, 0),
			extrcc = quantities::set_quantities(hk[grep("External Instrument", hkn), ], V, 0),
			lrefcc = quantities::set_quantities(hk[grep("Laser Reference", hkn), ], mW, 0),
			dpthcc = quantities::set_quantities(hk[grep("Depth", hkn), ], m, 0),
			tempcc = quantities::set_quantities(hk[grep("Temperature", hkn), ], `°C`, 0)
		)

		save(".LISSTi", file = system.file("R", "sysdata.rda", package = "lisst"))
		assign(".LISSTi", .LISSTi, envir = environment(lisst::lisst_reg))
		if(.LISSTi[[as.character(dat[1])]]$X) model <- paste0(model, "X")
		cat("\n")
		cat(paste0("LISST-", model), paste0("SN:", dat[1]), "successfully registered\n\n")
	}
	if(model == "200") {
		cat("To be added soon...")
	}
}

#'
#' For the LISST-100 series, extract only the relevant model information for the 
#' detector type.
#'

.getmodp <- function(linst) {
	lmodl      <- .LISSTm[[linst$mod]]
	lmodl$a0   <- lmodl$a0[linst$dty]
	lmodl$s0   <- lmodl$s0[[linst$dty]]
	lmodl$wang <- lmodl$wang[[linst$dty]]
	lmodl$binr <- lmodl$binr[[linst$dty]]
	lmodl
}

#' .LISSTm - A list of model dependent parameters list for each supported model.
#'
#' Variables:
#' pl    - Standard path length (m) of the instrument.
#' bnvar - Number of variables in the binary file.
#' anvar - Number of variables in the processed file.
#' nring - Number of detectors.
#' a0    - Minimum angle in air for the inner most ring (rad).
#' s0    - Minimum size range (µm) per inversion type (ss or rs) and per 
#' detector type (model dependent).
#' lvarn - Long variable names.
#' varun - Variable units.
#' aang  - In air angles (rad) for each detector ring (assinged to the correct 
#' bins), per detector type (model dependent).
#' wang  - In water angles (rad) for each detector ring(assinged to the correct 
#' bins), per detector type (model dependent).
#' binr  - Bin size range (µm) per inversion type (ss or rs) and per detector 
#' type (model dependent).

.LISSTm <- list(
	"100"  = list(
		mod   = "100",
		pl    = quantities::set_quantities(0.05, m, 0), 
		bnvar = 40, 
		anvar = 42, 
		nring = 32, 
		a0 = units::set_units(c(B = 0.1, C = 0.05) * pi / 180, rad),
		s0 = list(
			B = units::set_units(c(ss = 1.25, rs = 1.0), µm), 
			C = units::set_units(c(ss = 2.50, rs = 1.9), µm)
		),
		lvarn = c(paste("Bin", formatC(1:32, width = 2, flag = 0), sep = "_"),"Laser transmission","Battery voltage","External input 1","Laser reference","Depth","Temperature","Time1","Time2","Optical transmission","Beam attenuation"),
		varun = c(rep("ppm", 32), "mW", "V", "V", "mW", "m", "`°C`", "1", "1", "1", "1/m"),
		wang = list(
			B = units::set_units(matrix(c(
			rev(c(0.075,0.089,0.105,0.124,0.146,0.172,0.203,0.24,0.283,0.334,0.394,0.465,0.548,0.647,0.764,0.901,1.063,1.255,1.481,1.747,2.062,2.433,2.871,3.389,3.999,4.719,5.568,6.571,7.754,9.151,10.8,12.74)),
			rev(c(0.089,0.105,0.124,0.146,0.172,0.203,0.24,0.283,0.334,0.394,0.465,0.548,0.647,0.764,0.901,1.063,1.255,1.481,1.747,2.062,2.433,2.871,3.389,3.999,4.719,5.568,6.571,7.754,9.151,10.8,12.74,15.04)),
			rev(c(0.082,0.096,0.114,0.134,0.158,0.187,0.221,0.26,0.307,0.362,0.428,0.505,0.596,0.703,0.829,0.979,1.155,1.363,1.609,1.898,2.24,2.643,3.119,3.681,4.344,5.126,6.049,7.138,8.424,9.941,11.73,13.84))
			), ncol = 3, dimnames = list(paste("Bin", formatC(1:32, width = 2, flag = 0), sep = "_"), c("Min", "Max", "Med"))) * pi / 180, rad), 
			C = units::set_units(matrix(c(
			rev(c(0.038,0.044,0.052,0.062,0.073,0.086,0.102,0.12,0.141,0.167,0.197,0.232,0.274,0.324,0.382,0.451,0.532,0.627,0.74,0.874,1.031,1.217,1.436,1.694,1.999,2.359,2.784,3.286,3.877,4.575,5.399,6.371)),
			rev(c(0.044,0.052,0.062,0.073,0.086,0.102,0.12,0.141,0.167,0.197,0.232,0.274,0.324,0.382,0.451,0.532,0.627,0.74,0.874,1.031,1.217,1.436,1.694,1.999,2.359,2.784,3.286,3.877,4.575,5.399,6.371,7.519)),
			rev(c(0.041,0.048,0.057,0.067,0.079,0.093,0.11,0.13,0.154,0.181,0.214,0.252,0.298,0.351,0.415,0.489,0.578,0.682,0.804,0.949,1.12,1.322,1.56,1.841,2.172,2.563,3.025,3.569,4.212,4.97,5.865,6.921))
			), ncol = 3, dimnames = list(paste("Bin", formatC(1:32, width = 2, flag = 0), sep = "_"), c("Min", "Max", "Med"))) * pi / 180, rad)
		),
		binr = list(
			B = list(
				ss = units::set_units(matrix(c(
				c(1.25,1.48,1.74,2.05,2.42,2.86,3.38,3.98,4.7,5.55,6.55,7.72,9.12,10.8,12.7,15,17.7,20.9,24.6,29.1,34.3,40.5,47.7,56.3,66.5,78.4,92.6,109,129,152,180,212),
				c(1.48,1.74,2.05,2.42,2.86,3.38,3.98,4.7,5.55,6.55,7.72,9.12,10.8,12.7,15,17.7,20.9,24.6,29.1,34.3,40.5,47.7,56.3,66.5,78.4,92.6,109,129,152,180,212,250),
				c(1.36,1.6,1.89,2.23,2.63,3.11,3.67,4.33,5.11,6.03,7.11,8.39,9.9,11.7,13.8,16.3,19.2,22.7,26.7,31.6,37.2,43.9,51.9,61.2,72.2,85.2,101,119,140,165,195,230)
				), ncol = 3, dimnames = list(paste("Bin", formatC(1:32, width = 2, flag = 0), sep = "_"), c("Min", "Max", "Med"))), µm), 
				rs = units::set_units(matrix(c(
				c(1,1.18,1.39,1.64,1.94,2.29,2.7,3.19,3.76,4.44,5.24,6.18,7.29,8.61,10.2,12,14.1,16.7,19.7,23.2,27.4,32.4,38.2,45.1,53.2,62.8,74.1,87.4,103,122,144,169),
				c(1.18,1.39,1.64,1.94,2.29,2.7,3.19,3.76,4.44,5.24,6.18,7.29,8.61,10.2,12,14.1,16.7,19.7,23.2,27.4,32.4,38.2,45.1,53.2,62.8,74.1,87.4,103,122,144,169,200),
				c(1.09,1.28,1.51,1.79,2.11,2.49,2.93,3.46,4.09,4.82,5.69,6.71,7.92,9.35,11,13,15.4,18.1,21.4,25.2,29.8,35.2,41.5,49,57.8,68.2,80.5,94.9,112,132,156,184)
				), ncol = 3, dimnames = list(paste("Bin", formatC(1:32, width = 2, flag = 0), sep = "_"), c("Min", "Max", "Med"))), µm) 
			), 
			C = list(
				ss = units::set_units(matrix(c(
				c(2.5,2.95,3.48,4.11,4.85,5.72,6.75,7.97,9.4,11.1,13.1,15.4,18.2,21.5,25.4,30,35.4,41.7,49.2,58.1,68.6,80.9,95.5,113,133,157,185,218,258,304,359,424),
				c(2.95,3.48,4.11,4.85,5.72,6.75,7.97,9.4,11.1,13.1,15.4,18.2,21.5,25.4,30,35.4,41.7,49.2,58.1,68.6,80.9,95.5,113,133,157,185,218,258,304,359,424,500),
				c(2.72,3.2,3.78,4.46,5.27,6.21,7.33,8.65,10.2,12.1,14.2,16.8,19.8,23.4,27.6,32.5,38.4,45.3,53.5,63.1,74.5,87.9,104,122,144,170,201,237,280,331,390,460)
				), ncol = 3, dimnames = list(paste("Bin", formatC(1:32, width = 2, flag = 0), sep = "_"), c("Min", "Max", "Med"))), µm), 
				rs = units::set_units(matrix(c(
				c(1.9,2.25,2.65,3.13,3.69,4.35,5.14,6.06,7.15,8.44,9.96,11.8,13.9,16.4,19.3,22.8,26.9,31.8,37.5,44.2,52.2,61.6,72.7,85.7,101,119,141,166,196,232,273,322),
				c(2.25,2.65,3.13,3.69,4.35,5.14,6.06,7.15,8.44,9.96,11.8,13.9,16.4,19.3,22.8,26.9,31.8,37.5,44.2,52.2,61.6,72.7,85.7,101,119,141,166,196,232,273,322,381),
				c(2.07,2.44,2.88,3.4,4.01,4.73,5.58,6.59,7.77,9.17,10.8,12.8,15.1,17.8,21,24.8,29.2,34.5,40.7,48,56.7,66.9,78.9,93.1,110,130,153,181,213,252,297,350)
				), ncol = 3, dimnames = list(paste("Bin", formatC(1:32, width = 2, flag = 0), sep = "_"), c("Min", "Max", "Med"))), µm) 
			)
		)
	), 
	"200"  = list(
		mod   = "200",
		pl    = quantities::set_quantities(0.025, m, 0), 
		bnvar = 59, 
		anvar = 61, 
		nring = 36, 
		a0 = units::set_units(c(A = NA) * pi / 180, rad),
		s0 = list(
			A = units::set_units(c(ss = 1.00, rs = 1.00), µm)
		),
		lvarn = c(paste("Bin", formatC(1:36, width = 2, flag = 0), sep = "_"),"Laser transmission","Battery voltage","External input 1","Laser reference","Depth","Temperature","Year","Month","Day","Hour","Minute","Second","External input 2","Mean Diameter","Total Volume Concentration","Relative Humidity", "Accelerometer X", "Accelerometer Y", "Accelerometer Z", "Raw pressure 1","Raw pressure 2","Ambient Light", "NU", "Optical transmission","Beam attenuation"),
		varun = c(rep("ppm", 36), "mW", "V", "V", "mW", "m", "`°C`", "1", "1", "1", "hr", "min", "s", "V", "µm", "ppm", "`%`", "1", "1", "1", "1", "1", "1", "1", "1", "1/m"),
		wang = list(
			A = NULL
		),
		binr = list(
			A = list(
				ss = units::set_units(matrix(c(
				c(1,1.48,1.74,2.05,2.42,2.86,3.38,3.98,4.7,5.55,6.55,7.72,9.12,10.8,12.7,15,17.7,20.9,24.6,29.1,34.3,40.5,47.7,56.3,66.5,78.4,92.6,109,129,152,180,212,250,297,354,420),
				c(1.48,1.74,2.05,2.42,2.86,3.38,3.98,4.7,5.55,6.55,7.72,9.12,10.8,12.7,15,17.7,20.9,24.6,29.1,34.3,40.5,47.7,56.3,66.5,78.4,92.6,109,129,152,180,212,250,297,354,420,500),
				c(1.21,1.6,1.89,2.23,2.63,3.11,3.67,4.33,5.11,6.03,7.11,8.39,9.9,11.7,13.8,16.3,19.2,22.7,26.7,31.6,37.2,43.9,51.9,61.2,72.2,85.2,101,119,140,165,195,230,273,324,386,459)
				), ncol = 3, dimnames = list(paste("Bin", formatC(1:36, width = 2, flag = 0), sep = "_"), c("Min", "Max", "Med"))), µm),
				rs = units::set_units(matrix(c(
				c(1,1.48,1.74,2.05,2.42,2.86,3.38,3.98,4.7,5.55,6.55,7.72,9.12,10.8,12.7,15,17.7,20.9,24.6,29.1,34.3,40.5,47.7,56.3,66.5,78.4,92.6,109,129,152,180,212,250,297,354,420),
				c(1.48,1.74,2.05,2.42,2.86,3.38,3.98,4.7,5.55,6.55,7.72,9.12,10.8,12.7,15,17.7,20.9,24.6,29.1,34.3,40.5,47.7,56.3,66.5,78.4,92.6,109,129,152,180,212,250,297,354,420,500),
				c(1.21,1.6,1.89,2.23,2.63,3.11,3.67,4.33,5.11,6.03,7.11,8.39,9.9,11.7,13.8,16.3,19.2,22.7,26.7,31.6,37.2,43.9,51.9,61.2,72.2,85.2,101,119,140,165,195,230,273,324,386,459)
				), ncol = 3, dimnames = list(paste("Bin", formatC(1:36, width = 2, flag = 0), sep = "_"), c("Min", "Max", "Med"))), µm)
			)
		)
	)
)

# .LISSTi - A list of registered LISST instruments.
#
# Variables depend on model.
#
# On creation only:
#
# .LISSTi <- list()
#
# save(.LISSTm, .LISSTi, file = system.file("R", "sysdata.rda", package = "lisst"))

