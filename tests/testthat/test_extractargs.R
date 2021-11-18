#context("matchcall of diff_analysis")
#data(kostic2012crc)
#set.seed(1024)
#targetgroup <- "DIAGNOSIS"
#diffres <- diff_analysis(kostic2012crc, classgroup=targetgroup,
#		        mlfun="lda", filtermod="fdr",
#		        firstcomfun = "kruskal.test", 
#		        firstalpha=0.05, strictmod=TRUE,
#	                submin=3, subclwilc=TRUE,
#	                secondcomfun = "wilcox.test",
#	                secondalpha=0.01, lda=3)
#
#test_that("extract_args",{
#    expect_true(is(diffres@someparams, "list"))
#})
#
#context("extract_args function")
#test_that("extract_args",{
#    expect_equal("kruskal.test", MicrobiotaProcess:::extract_args(diffres, "firstcomfun"))
#    expect_equal("wilcox.test", MicrobiotaProcess:::extract_args(diffres, "secondcomfun"))
#    expect_equal("lda", MicrobiotaProcess:::extract_args(diffres, "mlfun"))
#    expect_equal(targetgroup, MicrobiotaProcess:::extract_args(diffres, "classgroup"))
#})
#
