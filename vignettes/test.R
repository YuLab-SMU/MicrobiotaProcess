source("../R/CountOrRatios.R")
otuda <- read.table("../data/otu_tax_table.txt",
                  sep="\t", header=T, row.names=1, check.names=F, comment.char="",skip=1)
tax <- otuda[,"taxonomy",drop=FALSE]
otuda$taxonomy <- NULL
otuda <- data.frame(t(otuda), check.names=F)
head(otuda[,1:7])
sample <- read.table("../data/sample_info.txt", header=T, row.names=1)
dat <- CountOrRatios(otuda, sample, countmode=TRUE,
                          percentmode=FALSE,
                          multiplenum=1,
                          rownamekeep=FALSE)
head(dat, 3)

