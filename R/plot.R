
# if xu is not defined the default unit will be used.

plot.lisst <- function(lo, xu, type, ...) {

	linst <- attr(lo, "linst")
	lmodl <- attr(lo, "lmodl")
	typ   <- attr(lo, "type") 
	if(missing(type))
		type <- attr(lo, "type")

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
		x <- lmodl$binr[[linst$ity]]
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
