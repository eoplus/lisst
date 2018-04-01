
# if xu is not defined the default unit will be used.

plot.lisst <- function(lo, xu, type, ...) {

	linst <- attr(lo, "linst")
	lmodl <- attr(lo, "lmodl")

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
			if(length(id) > 0) cp <- cp[-id, ]
			y <- y / cp
			plot.vsf(x, y, pf = TRUE, ...)
		} else {
			plot.vsf(x, y, pf = FALSE, ...)
		}
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
