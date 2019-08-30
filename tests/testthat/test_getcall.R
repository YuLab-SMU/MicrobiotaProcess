context("matchcall of diffAnalysis")
data(kostic2012crc)
set.seed(1024)
diffres <- diffAnalysis(kostic2012crc, class="DIAGNOSIS",
						mlfun="lda", filtermod="fdr",
						firstcomfun = "kruskal.test", 
						firstalpha=0.05, strictmod=TRUE,
						submin=3, subclwilc=TRUE,
						secondcomfun = "wilcox.test",
						secondalpha=0.01, lda=3)

test_that("getcall",{
    expect_true(is(diffres@call, "call"))
})

context("getcall function")
test_that("getcall",{
   	expect_equal("kruskal.test", MicrobiotaProcess:::getcall(diffres, "firstcomfun"))
    expect_equal("wilcox.test", MicrobiotaProcess:::getcall(diffres, "secondcomfun"))
    expect_equal("lda", MicrobiotaProcess:::getcall(diffres, "mlfun"))
})

