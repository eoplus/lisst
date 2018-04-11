## ----echo=FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE)

## ------------------------------------------------------------------------
# Register the LISST instrument:
path  <- system.file("extdata", package = "lisst")
model <- 100
lisst_reg(model, path)

# Read a LISST data file:
sn    <- 1298
pl    <- 0.05
yr    <- 2018 
zscat <- system.file("extdata", "bg_20180326.asc", package = "lisst")
flp   <- system.file("extdata", "DN_27_rs.asc", package = "lisst") # LISST-SOP processed file
flb   <- system.file("extdata", "DN_27.DAT", package = "lisst")    # LISST binary file

lop   <- read_lisst(flp, sn, pl, zscat, yr)
lob   <- read_lisst(flb, sn, pl, zscat, yr)

# Subset by sample, deth or time:
lob_ss <- lob[1:3, ]         # First three samples
lob_ss <- lob['0|5', ]       # First five meters
lob_ss <- lob['2018-03', ]   # Samples in march 2018

# Conversion between data types:
lop_pnc <- lget(lop, 'pnc')
lob_cal <- lget(lob, 'cal')

# Concatenation:
lop_d <- c(lop, lop_pnc)

# Statistics:
#lob_mn <- lstat(lob, mean)
#lob_md <- lstat(lob, median)
#breaks <- c(0, 5, 10, 15)
#lob_mn <- lgetstat(lob, mean, s = breaks) # Breaks by sample index
#lob_mn <- lgetstat(lob, mean, d = breaks) # Breaks by sample depths
#breaks <- c("2018-01", "2018-02", "2018-03", "2018-04")
#lob_mn <- lgetstat(lob, mean, t = breaks) # Breaks by sample time
 
# Plots:
# par(mfcol = c(2, 2))
# plot(lob[40:144, ])
# plot(lop[40:144, ], type = 'pnc')
# plot(lob[40:144, ], type = 'pf', xu = 'degree')
# plot(lop[40:144, ], type = 'vol')

## ------------------------------------------------------------------------
#path  <- system.file("extdata", package = "lisst")
#model <- 100
#lisst_reg(model, path)

## ------------------------------------------------------------------------
#sn    <- 1298
#pl    <- 0.05
#yr    <- 2018 
#zscat <- system.file("extdata", "bg_20180326.asc", package = "lisst")

# For a processed file:
#flp   <- system.file("extdata", "DN_27_rs.asc", package = "lisst")
#lop   <- read_lisst(flp, sn, pl, zscat, yr)

# For a binary file:
#flb   <- system.file("extdata", "DN_27.DAT", package = "lisst")
#lob   <- read_lisst(flb, sn, pl, zscat, yr)

## ------------------------------------------------------------------------
#lom <- read_lisst(flp, pl = pl, yr = yr, model = "100CX")

## ------------------------------------------------------------------------
#lop[15, , drop = FALSE]
#lob[15, , drop = FALSE]

## ------------------------------------------------------------------------
# Converting from vol conc to n conc:
#lop[15, 1:3, drop = FALSE]
#lgetpnc(lop)[15, 1:3, drop = FALSE]

# Converting from cal to VSF:
#lob[15, 1:3, drop = FALSE]
#lgetvsf(lob)[15, 1:3, drop = FALSE]

## ------------------------------------------------------------------------
# Extract the Junge slope for a PSD:
#lgetfit(lop)[15:17]

## ------------------------------------------------------------------------
#lov <- lgetvsf(lob)
#par(mfcol = c(2, 2))
#plot(lov[40:144, ])
#plot(lop[40:144, ], type = 'pnc')
#plot(lov[40:144, ], type = 'pf', xu = 'degree')
#plot(lop[40:144, ], type = 'vol')

