## ----echo=FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE)

## ------------------------------------------------------------------------
path  <- system.file("extdata", package = "lisst")
model <- 100
lisst_reg(model, path)

## ------------------------------------------------------------------------
sn    <- 1298
pl    <- 0.05
yr    <- 2018 
zscat <- system.file("extdata", "bg_20180326.asc", package = "lisst")

# For a processed file:
flp   <- system.file("extdata", "DN_27.asc", package = "lisst")
lop   <- read_lisst(flp, sn, pl, zscat, yr)

# For a binary file:
flb   <- system.file("extdata", "DN_27.DAT", package = "lisst")
lob   <- read_lisst(flb, sn, pl, zscat, yr)

## ------------------------------------------------------------------------
lom <- read_lisst(flp, pl = pl, yr = yr, model = "100CX")

## ------------------------------------------------------------------------
lop[15, , drop = FALSE]
lob[15, , drop = FALSE]

## ------------------------------------------------------------------------
# Converting from vol conc to n conc:
lop[15, 1:3, drop = FALSE]
lgetpnc(lop)[15, 1:3, drop = FALSE]

# Converting from cal to VSF:
lob[15, 1:3, drop = FALSE]
lgetvsf(lob)[15, 1:3, drop = FALSE]

## ------------------------------------------------------------------------
lov <- lgetvsf(lob)
par(mfcol = c(2, 2))
plot(lov[40:144, ])
plot(lov[40:144, ], type = "pf")
plot(lov[40:144, ], xu = "degree")
plot(lov[40:144, ], xu = "degree", type = "pf")

