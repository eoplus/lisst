
devtools::load_all()

model = "100"
path = "inst/extdata"
lisst_reg(model, path)

fl    <- "inst/extdata/DN_27.DAT"
zscat <- "inst/extdata/bg_20180326.asc"
sn    <- "1298"
pl    <- 0.05
yr    <- 2018

fl    <- "inst/extdata/DN_27.DAT"
x     <- read_lisst(fl, sn, pl, zscat, yr, out = 'vsf')
donkmeer_bin <- x[-c(1:13), ]
rownames(donkmeer_bin) <- 1:nrow(donkmeer_bin)
devtools::use_data(donkmeer_bin, overwrite = TRUE)

fl    <- "inst/extdata/DN_27_rs.asc"
x     <- read_lisst(fl, sn, pl, zscat, yr, out = 'vol')
donkmeer_pro <- x[-c(1:13), ]
rownames(donkmeer_pro) <- 1:nrow(donkmeer_pro)
devtools::use_data(donkmeer_pro, overwrite = TRUE)

###################################################
# Fictional data for testing and exemple purposes #
###################################################

# LISST-100X:

base <- c(1, 1.5, 2, 2.5, 3, 4, 5, 5.5, 5, 5, 4, 3, 2.5, 2, 1.1, 1, 0.9, 1.1, 
	1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 7.1, 6.9, 6, 4, 2)
base    <- base / sum(base)
depth_p <- seq(0, 1, 0.1)^2
conc    <- 100 
depth   <- seq(0, 1, 0.05)
depth   <- c(depth, rev(depth)[-1])
depth   <- depth + c(0.032, rnorm(length(depth)-2, mean = 0, sd = 0.0125), 0.041)
depth[depth > 1] <- 1
fact    <- approx(x = seq(0, 1, 0.1), y = depth_p, xout = depth)$y * conc
ring    <- matrix(base, ncol = 32, nrow = length(depth), byrow = T)
for(i in 1:nrow(ring)) ring[i, ] <- ring[i, ] * fact[i]
OptTrans <- (depth^-0.3) / 3
RLaser   <- rep(1.75178, length(depth))
TLaser   <- RLaser * OptTrans
Battery  <- seq(7.73, 7.33, by = -0.01)
ExtI1    <- rep(0.01002, length(depth))
Temperature <- 10.2 + rnorm(41, 0, 0.1)
BeamAtt  <- -log(OptTrans)/0.05
Depth 	 <- depth * 10 
base <- cbind(ring, TLaser, Battery, ExtI1, RLaser, Depth, Temperature, OptTrans, BeamAtt)
base <- as.data.frame(base)
lmodl <- lisst:::.getmodp(list(mod = "100", dty = "C"))
linst <- list(X = TRUE)
linst$mod <- '100C'
zscatd <- rep(as.numeric(NA), lmodl$bnvar)
zscatd <- as.data.frame(matrix(zscatd, nrow = 1))
names(zscatd) <- lmodl$lvarn[1:lmodl$bnvar]
colnames(base) <- lmodl$lvarn[-c(39:40)]
id <- c(1:38)
for(i in id) {
	units(base[, i]) <- lmodl$varun[i]
	if(i < lmodl$bnvar) units(zscatd[, i]) <- 1
}
units(base[, 39]) <- lmodl$varun[41]
units(base[, 40]) <- lmodl$varun[42]
ti <- seq(as.POSIXct("2018-04-03 10:05:37", tz = 'UTC'), as.POSIXct("2018-04-03 10:10:37", tz = 'UTC'), length.out = 41)
rownames(base) <- format(ti, "%Y-%m-%d %H:%M:%OS1 %Z")
fprofile <- structure(base, type = 'vol', lproc = list(ity = 'rs'), 
		linst = linst, lmodl = lmodl, zscat = zscatd, 
		class = c("lisst", "data.frame"))
devtools::use_data(fprofile, overwrite = TRUE)



