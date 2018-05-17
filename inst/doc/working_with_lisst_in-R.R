## ----echo=FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE)
library(units)
units_options(negative_power = TRUE, parse = TRUE, group = c("(", ")"))

## ------------------------------------------------------------------------
library(lisst)

fl200 <- system.file("extdata", "sp_april_rs.csv", package = "lisst")
l200p <- read_lisst(fl200)

## ------------------------------------------------------------------------
head(l200p, n = 2)

## ---- fig.height = 7, fig.width = 6, fig.align="center"------------------
lhov(l200p, by = 'time', norm = F)

## ---- fig.height = 7, fig.width = 6, fig.align="center"------------------
l200p <- l200p['2018-04-19 07:30/2018-04-19 14:00', ]
lhov(l200p, by = 'time', norm = F)

## ---- fig.height = 4, fig.width = 5, fig.align="center"------------------
plotBins(l200p, bins = list(1:14, 25:36), by = 'time', col = c("blue", "red"))
binr  <- lbinr(l200p)
r1    <- range(binr[1:14,])
r2    <- range(binr[25:36,])
ltext <- c(eval(substitute(expression(x-y~mu*m), list(x = r1[1], y = r1[2]))),
           eval(substitute(expression(x-y~mu*m), list(x = r2[1], y = r2[2]))))
legend('topleft', legend = ltext, lty = 1, col = c("blue", 'red'), bty = 'n', cex = 0.8)

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
#lget(lop, 'pnc')[15, 1:3, drop = FALSE]

# Converting from cal to VSF:
#lob[15, 1:3, drop = FALSE]
#lget(lob, 'vsf')[15, 1:3, drop = FALSE]

## ------------------------------------------------------------------------
# Extract the Junge slope for a PSD:
#lfit(lop)[15:17]

## ------------------------------------------------------------------------
#lov <- lgetvsf(lob)
#par(mfcol = c(2, 2))
#plot(lov[40:144, ])
#plot(lop[40:144, ], type = 'pnc')
#plot(lov[40:144, ], type = 'pf', xu = 'degree')
#plot(lop[40:144, ], type = 'vol')

## ---- fig.height = 7, fig.width = 7--------------------------------------
#lov <- lget(lob, 'vsf')
#lhov(lob, by = 'sample')
#lhov(lob, by = 'sample', norm = FALSE)

