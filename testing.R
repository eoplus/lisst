
library(units)
units_options(negative_power = TRUE, parse = TRUE, group = c("(", ")"))
source("R/register.R")
source("R/import.R")
source("R/manipulate.R")
source("R/process.R")
source("R/plot.R")

model = "100"
path = "/home/alexandre/Documents/research/belgian_campaigns/field_data/data/lisst100/calfiles_sn_1298"

lisst_reg(model, path)

fl    <- "test/DN_27.DAT"
zscat <- "test/bg_20180326.asc"
sn <- 1298

testBin1 <- lisst_bin(fl, zscat, sn, out = "raw")
testBin2 <- lisst_bin(fl, zscat, sn, out = "cor")
testBin3 <- lisst_bin(fl, zscat, sn, out = "cal")
testvsf1 <- getvsf(testBin3)

plot(testvsf)

fl    <- "test/SP_37.DAT"
zscat <- "test/bg_20180307.asc"
lbin  <- lisst_bin(fl, zscat, sn, out = "cal")
testvsf2 <- getvsf(lbin)

testvsf <- rbind(testvsf1[20:144], testvsf2[15:45,][-13, ])

plot(testvsf)



fl <- "test/DN_27.asc"
zscat <- "test/bg_20180326.asc"
sn <- 1298

testAsc <- lisst_asc(fl, sn)



fl    <- "test/DN_27.log"
zscat <- "test/bg_20180326.asc"
sn <- 1298
testAsc <- lisst_asc(fl, sn)
