context("mp_cal_venn and getvennlist work")
test_that("getvennlist work",{
	dir <- system.file("extdata", package="MicrobiotaProcess")
	otuda <- read.table(file.path(dir,"otu_tax_table.txt"),
						row.names=1, check.names=FALSE, skip=1,
						comment.char="", header=TRUE, sep="\t")
	otuda <- otuda[vapply(otuda,is.numeric,logical(1))]
	otuda <- data.frame(t(otuda))
	sampleda <- read.table(file.path(dir, "sample_info.txt"), 
						   header=TRUE, row.names=1)
	vennlist1 <- get_vennlist(obj=otuda, sampleinfo=sampleda, factorNames="group")
	testph <- mp_import_qiime(file.path(dir,"otu_tax_table.txt"),
						   file.path(dir, "sample_info.txt"))
    vennlist2 <- testph %>% 
                 mp_cal_venn(.abundance=Abundance, .group=group, action="get", force=TRUE) %>% 
                 dplyr::pull(vennOfgroup, name=group)
	expect_equal(vennlist1, vennlist2)
})
