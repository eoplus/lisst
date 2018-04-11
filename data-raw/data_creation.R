

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
devtools::use_data(donkmeer_bin, overwrite = TRUE)

fl    <- "inst/extdata/DN_27_rs.asc"
x     <- read_lisst(fl, sn, pl, zscat, yr, out = 'vol')
donkmeer_pro <- x[-c(1:13), ]
devtools::use_data(donkmeer_pro, overwrite = TRUE)
