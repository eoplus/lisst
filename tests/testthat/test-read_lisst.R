context("Import LISST data")

test_that("read LISST-100X binary file works", {
	flSOP    <- system.file("extdata", "DN_27.log", package = "lisst")
	binSOP   <- read.table(flSOP, header = F)
	for(i in 1:38) units(binSOP[, i]) <- 1
	fltest   <- system.file("extdata", "DN_27.DAT", package = "lisst")
	binlisst <- read_lisst(fltest, sn = '1298', pl = 0.05, yr = 2018, out = 'raw')
  for(i in 1:40) expect_equal(binlisst[, i], binSOP[, i])
})

test_that("read LISST-200X binary file works", {
	flSOP    <- system.file("extdata", "sp_april.rtx", package = "lisst")
	binSOP   <- read.table(flSOP, header = F, sep = ',')
	binSOP[binSOP > 40950] <- binSOP[binSOP > 40950] - 65536
	binSOP[, 1:36] <- binSOP[, 1:36] / 10
	for(i in c(1:42, 49:59)) units(binSOP[, i]) <- 1
	fltest   <- system.file("extdata", "sp_april.RBN", package = "lisst")
	binlisst <- read_lisst(fltest, sn = '2028', pl = 0.025, out = 'raw')
  for(i in 1:59) expect_equal(binlisst[, i], binSOP[, i])
})
