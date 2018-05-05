
# Pure water optical properties at 670 nm:
# 
# References:
# Water Optical Properties Processor (WOPP)

aw670 <- 0.439       # 1/m
bw670 <- 5.808404e-4 # 1/m


# Capitalization function
# From ?toupper examples

.cap <- function(x) {
	s <- strsplit(x, " ")[[1]]
	paste(toupper(substring(s, 1, 1)), substring(s, 2),
		sep = "", collapse = " ")
}


