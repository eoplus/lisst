
#' Auxiliary plotting functions

.xaxn <- function(x, by) {
	switch(by,
	sample = list("Sample Number", 1:nrow(x)),
	depth  = list(units::make_unit_label("'Depth'", x$Depth), x$Depth),
	time   = list("Time", x$Time))
}

.yaxn <- function(x) {
	lmodl <- attr(x, "lmodl")
	ity   <- attr(x, "lproc")$ity
	switch(attr(x, "type"),
		raw = ,
		cor = ,
		cal = list("'Bins'", set_units(1:lmodl$nring, 1)),
		vsf = list("'Scattering Angle'", lmodl$wang[, 3]),
		vol = ,
		pnc = list("'Size'", lmodl$binr[[ity]][, 3])
	)
}

.zaxn <- function(x) {
	switch(attr(x, "type"),
	raw = ,
	cor = "'Digital Number'",
	cal = "'Laser Power'",
	vsf = "'Volume Scattering Function'",
	vol = "'Volume Concentration'",
	pnc = "'Number Concentration'")
}

.map2color <- function(x, pal, limits){
    if(missing(limits)) limits <- range(x, na.rm = T)
    pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal)+1), all.inside = TRUE)]
}


#' Hovmöller diagram for lisst objects
#'
#' Create Hovmöller diagram for lisst objects based on sample number, 
#' measurement time or depth.
#'
#' @param x      A lisst object.
#' @param by     Ordinate axis dimension. Must be one of 'sample', 'time' or 
#'               'depth'.
#' @param nbins  If the lisst object is not regular in the dimension specified 
#'               parameter by, nbins controls the binning of the data to achieve 
#'               a regular series appropriate for raster representation.
#' @param norm   Logic. Should the magnitude values de normalized? See details.
#' @param legend Logic. Should a legend be added to the plot?
#' @param col    A color pallete to be used. If not specified a default pallete 
#'               will be used.
#' @param yu     The units for the abscissa. See details.
#' @param zu     The units for the magnitude values. See details.
#'
#' @details
#' If the lisst object is not regularly spaced in the chosen dimension, it will 
#' be binned for raster representation. The nbins parameter will control the 
#' number of bins, with all bins having the same interval size. The sample 
#' dimension will always be regular.
#'
#' Normalization is performed by the sum of all detector values. It is provided 
#' to aid visual separation from changes in concentration and distribution. If 
#' TRUE a top pannel is added to the plot, with a line plot indicating the 
#' magnitude of the sum of detector values.
#'
#' If legend is set to TRUE a bottom pannel will be added to the plot with a 
#' magnitude scale. The axis label and units are affected by the zlab and zu 
#' parameters.
#'
#' Note that unit conversion is performed 'automatcaly' using the units package. 
#' If there is no standard conversion between the default and the request units
#' the function will exit with an error. See \code{link{units}} for details.
#'
#' @examples
#' flp <- system.file("extdata", "DN_27_rs.asc", package = "lisst")
#' model <- "100CX" 
#' lop <- read_lisst(flp, model = model)
#' lhov(lop, by = 'sample')
#' lhov(lop, by = 'sample', norm = FALSE)
#' lhov(lop, by = 'sample', norm = FALSE, legend = FALSE)
#'
#' @export

by = 'sample'
nbins = "pretty"
norm = TRUE
legend = TRUE


lhov <- function(x, by = 'sample', nbins = "pretty", norm = TRUE, legend = TRUE, col, xlab, ylab, 
	zlab, yu, zu) {

	lmodl  <- attr(x, "lmodl")
	lty    <- attr(x, "type")
	xm     <- do.call(rbind, x[, 1:lmodl$nring])
	if(lty == 'vol' | lty == 'pnc') {
		.na <- set_errors(as.numeric(NA), 0)
		units(.na) <- units(xm)
		xm[drop_quantities(xm) == 0] <- .na
	}
	if(legend && norm) {
		layout(matrix(c(1, 2, 3), ncol = 1), heights = c(1, 4, 1))
	} else if(legend) {
		layout(matrix(c(1, 2), ncol = 1), heights = c(4, 1))
	} else if(norm) {
		layout(matrix(c(1, 2), ncol = 1), heights = c(1, 4))
	}

	xaxn <- .xaxn(x, by)
	if(missing(xlab)) xlab <- xaxn[[1]]
	yaxn <- .yaxn(x)
	if(!missing(yu)) units(yaxn[[2]]) <- yu
	if(missing(ylab)) ylab <- units::make_unit_label(yaxn[[1]], yaxn[[2]])
	zaxn <- .zaxn(x)
	if(!missing(zu)) xm <- set_units(xm, zu)
	if(missing(zlab)) {
		if(norm) zlab <- paste("Normalized", parse(text = zaxn))
		else zlab <- units::make_unit_label(zaxn, x[, 1])
	}

	if(!all(diff(xaxn[[2]]))) stop('Binning for irregular dimension not implemented...')

	if(missing(col)) {
		col <- colorRampPalette(c("#011F4B", "#03396C", "#005B96", "#6497B1", "#B3CDE0"))(256)
	}

	if(norm) {
		xs <- xm[1, ]
		for(i in 1:ncol(xm)) xs[i] <- sum(xm[, i])
		xm   <- xm / rep(xs, each = nrow(xm))
		par(mar = c(0, 5, 0, 1))
		nylab <- units::make_unit_label('Total', x[, 1])
		plot(xaxn[[2]], drop_quantities(xs), axes = F, xlab = "", type = "l", xaxs = 'i', 
			ylab = nylab, log = 'y')
		axis(2)
	}

	par(mar = c(5, 5, 2, 1))
	plot(NA, xlim = c(1, ncol(xm) + 1), ylim = c(1, lmodl$nring + 1), bty = "l", xaxs = 'i', 
		yaxs = 'i', xaxt = 'n', yaxt = 'n', xlab = xlab, ylab = ylab)
	id <- c(seq(1, lmodl$nring, by = 4), lmodl$nring)
	axis(2, at = id+0.5, labels = round(yaxn[[2]][id], 4))
	axis(1, at = axTicks(1) + 0.5, labels = xaxn[[2]][axTicks(1)])
	xm   <- log10(drop_quantities(xm))
	xrst <- matrix(.map2color(xm, col), ncol = ncol(xm))[nrow(xm):1, ]
	rasterImage(xrst, xleft = 1, ybottom = 1, xright = ncol(xm) + 1, ytop = lmodl$nring + 1, 
		interpolate = F)

	if(legend) {
		par(mar = c(5, 8, 1, 4))
		plot(NA, xlim = range(xm, na.rm = T), ylim = c(0, 1), xaxs = 'i', yaxs = 'i', xaxt = 'n', 
			yaxt = 'n', xlab = zlab, ylab = "")
		digits <- as.numeric(sapply(strsplit(formatC(10^axTicks(1), format = 'e'), "e"), 
			'[[', 2))
		digits[digits > 0] <- 0
		axis(1, at = axTicks(1), labels = round(10^axTicks(1), abs(digits)))
		xrst <- t(as.raster(.map2color(seq(min(xm, na.rm = T), max(xm, na.rm = T), length.out = length(col)), col)))
		rasterImage(xrst, xleft = min(xm, na.rm = T), ybottom = 0, xright = max(xm, na.rm = T), ytop = 1, interpolate = T)
	}


# If not per sample, aggregate to create a regular serires.
#	if(by != 'sample') {
#		by <- x[, .cap(by)]
# 		For pretty bins...
# 		range(by) / length(x)
#		seq(by[1], by[length(by)], length.out = nbins)
#		if less than 5 days, if less than 5 hours, per 5 min, if less tha 5 min, per 5 seconds.
#	}
}


# if xu is not defined the default unit will be used.

plot.lisst <- function(lo, xu, type = "rings", ...) {

	linst <- attr(lo, "linst")
	lmodl <- attr(lo, "lmodl")
	typ   <- attr(lo, "type")

	if(missing(type))
		type <- attr(lo, "type")

	if(type == 'rings') {
		
	}


	if(type == "vsf" || type == "pf") {
		if(attr(lo, "type") != "vsf")
			stop("vsf or pf plot requires a vsf type lisst object", call. = FALSE)
		if(!missing(xu) && !(xu == "degree" || xu == "radian"))
			stop("xu for vsf plots must be 'degree' or 'rad'", call. = FALSE)
		x <- lmodl$wang[, 3]
		if(!missing(xu)) units(x) <- xu
		id <- which(lo[, 1] == set_units(0, 1/m/sr))
		y <- set_units(as.matrix(lo[, 1:lmodl$nring]), 1/m/sr)
		if(length(id) > 0) y <- y[-id, ]
		if(type == "pf") {
			cp <- lo[, "Beam attenuation"]
			aw670 <- units::set_units(0.439, 1/m)
			bw670 <- units::set_units(0.0005808404, 1/m)
			cp <- cp - ((aw670 + bw670) * as.numeric(lmodl$pl))
			if(length(id) > 0) cp <- cp[-id, ]
			y <- y / cp
			plot.vsf(x, y, pf = TRUE, ...)
		} else {
			plot.vsf(x, y, pf = FALSE, ...)
		}
	}
	if(type == 'vol' || type == 'pnc' || type == 'psd') {
		x <- lmodl$binr[[attr(lo, "lproc")]]
		if(!missing(xu)) units(x) <- xu
		id <- which(as.numeric(lo[, 1]) == 0)
		if((type == 'vol' && typ == 'pnc') || type == 'psd') lo <- lgetvol(lo)
		if(type == 'pnc' && typ == 'vol') lo <- lgetpnc(lo)
		y <- as.matrix(lo[, 1:lmodl$nring, drop = FALSE])
		units(y) <- units(lo[, 1])
		y2 <- NULL
		if(type == 'psd' || type == 'pnc') {
			lo <- lgetpnc(lo)
			y2 <- as.matrix(lo[, 1:lmodl$nring, drop = FALSE])
			units(y2) <- units(lo[, 1])
		}
		if(length(id) > 0) {
			y <- y[-id, ]
			if(!is.null(y2)) y2 <- y2[-id, ]
		}
		plot.psd(x, y, y2, type, ...)

	}
}

plot.psd <- function(x, y, y2, type, ...) {
	par(mar = c(5, 5, 4, 5))
	yn <- "'Particle Volume Concentration'"
	ylab <- units::make_unit_label(yn, y)
	if(type == 'pnc') {
		yn <- "'Particle Number Concentration'"
		ylab <- units::make_unit_label(yn, y2)
	}
	xlab <- units::make_unit_label("'Particle size'", x)
	ylim <- range(y, na.rm = T)
	if(type == 'pnc') ylim <- range(y2, na.rm = T)
	plot(NA, xlim = range(c(x[1, 1], x[, 2])), ylim = ylim, log = "xy", 
	xlab = xlab, ylab = ylab, bty = "l", yaxs = "i", ...)
	if(type == 'psd' || type == 'vol') {
		for(i in 1:ncol(y)) {
			rect(xleft = x[i, 1], ybottom = ylim[1], xright = x[i, 2], 
			ytop = y[1, i])
		}
	} else {
		for(i in 1:nrow(y2)) lines(x[, 3], y2[i, ], lty = 1)
	}
	if(type == 'psd') {
		par(new = T)
		plot(NA, xlim = range(c(x[1, 1], x[, 2])), ylim = range(y2), 
		log = "xy", xlab = "", ylab = "", bty = "n", main = "", xaxt = "n", 
		yaxt = "n")
		bor <- par("usr")
		lines(c(10^bor[2], 10^bor[2]), c(10^bor[3], 10^bor[4]), col = "red")
		axis(4, col = "red", col.axis = "red")
		lines(x[, 3], y2[1, ], lty = 1, col = "red")
		yn <- "'Particle Number Concentration'"
		ylab <- units::make_unit_label(yn, y2)
		mtext(ylab, side = 4, line = 3, col = "red")
	}
}

plot.vsf <- function(x, y, pf = FALSE, ...) {
	par(mar = c(5, 5, 3, 2))
	ylim <- range(y, na.rm = T)
	yn <- "'Particle Volume Scattering Function'"
	if(pf) yn <- "'Particle Phase Function'"
	ylab <- units::make_unit_label(yn, y)
	xlab <- units::make_unit_label("'In water scattering angle'", x) 
	plot(x, y[1, ], type = "l", log = "xy", ylim = ylim, xlab = xlab, ylab = ylab, ...)
	if(nrow(y) > 1)
		for(i in 1:nrow(y))
			lines(x, y[i, ])
}
