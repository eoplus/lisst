context("lisst attributes")

test_that("retrieval of in water angles works", {
	expect_identical(lwang('100B'), lisst:::.LISSTm[['100']]$wang$B)
	expect_identical(lwang('100C'), lisst:::.LISSTm[['100']]$wang$C)
	expect_identical(lwang('200'), lisst:::.LISSTm[['200']]$wang$A)

	fltest <- system.file("extdata", "DN_27_rs.asc", package = "lisst")
	lob <- read_lisst(fltest, sn = '1298', pl = 0.05, yr = 2018, out = 'vol')
	expect_identical(lwang(lob), lisst:::.LISSTm[['100']]$wang$C)

	fltest <- system.file("extdata", "sp_april.RBN", package = "lisst")
	lob <- read_lisst(fltest, pl = 0.025, out = 'raw')
	expect_identical(lwang(lob), lisst:::.LISSTm[['200']]$wang$A)
})

test_that("retrieval of size ranges works", {
	expect_identical(lbinr('100B', 'ss'), lisst:::.LISSTm[['100']]$binr$B$ss)
	expect_identical(lbinr('100B', 'rs'), lisst:::.LISSTm[['100']]$binr$B$rs)
	expect_identical(lbinr('100C', 'ss'), lisst:::.LISSTm[['100']]$binr$C$ss)
	expect_identical(lbinr('100C', 'rs'), lisst:::.LISSTm[['100']]$binr$C$rs)
	expect_identical(lbinr('200X', 'ss'), lisst:::.LISSTm[['200']]$binr$A$ss)
	expect_identical(lbinr('200X', 'rs'), lisst:::.LISSTm[['200']]$binr$A$rs)
	expect_identical(lbinr('200X'), lisst:::.LISSTm[['200']]$binr$A$ss)

	fltest <- system.file("extdata", "DN_27_rs.asc", package = "lisst")
	lop <- read_lisst(fltest, sn = '1298', pl = 0.05, yr = 2018, out = 'vol')
	expect_identical(lbinr(lop), lisst:::.LISSTm[['100']]$binr$C$rs)
	expect_identical(lbinr(lop, 'ss'), lisst:::.LISSTm[['100']]$binr$C$ss)

	fltest <- system.file("extdata", "DN_27.DAT", package = "lisst")
	lob <- read_lisst(fltest, sn = '1298', pl = 0.05, yr = 2018, out = 'raw')
	expect_identical(lbinr(lob, 'ss'), lisst:::.LISSTm[['100']]$binr$C$ss)

	fltest <- system.file("extdata", "sp_april.RBN", package = "lisst")
	lob <- read_lisst(fltest, pl = 0.025, out = 'raw')
	expect_identical(lbinr(lob), lisst:::.LISSTm[['200']]$binr$A$ss)
})

test_that("lbinr catch reasonable errors", {
	expect_error(lbinr('000'), "model must be one of '100(X)B', '100(X)C' or '200X'", fixed = T)
	expect_error(lbinr(100), "x must be a lisst or a character object")
	expect_error(lbinr('100XC'), "ity is not defined in the lisst object or needs more specification than just model")
	expect_error(lbinr('100XC', 'pp'), "ity must be either 'ss' or 'rs'")
})

test_that("lwang catch reasonable errors", {
	expect_error(lwang('000'), "model must be one of '100(X)B', '100(X)C' or '200X'", fixed = T)
	expect_error(lwang(100), "x must be a lisst or a character object")
})

