context("getvennlist work")
test_that("getvennlist work",{
	library(phyloseq)
	dir <- system.file("extdata", package="MicrobiotaProcess")
	otuda <- read.table(file.path(dir,"otu_tax_table.txt"),
						row.names=1, check.names=FALSE, skip=1,
						comment.char="", header=TRUE, sep="\t")
	otuda <- otuda[vapply(otuda,is.numeric,logical(1))]
	otuda <- data.frame(t(otuda))
	sampleda <- read.table(file.path(dir, "sample_info.txt"), 
						   header=TRUE, row.names=1)
	vennlist1 <- getvennlist(obj=otuda, sampleinfo=sampleda, 
							 factorNames="group")
	testph <- import_qiime(file.path(dir,"otu_tax_table.txt"),
						   file.path(dir, "sample_info.txt")) 
	vennlist2 <- getvennlist(testph,factorNames="group")
	expect_equal(vennlist1, vennlist2)
})
