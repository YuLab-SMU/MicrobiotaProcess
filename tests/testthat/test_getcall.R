context("matchcall of diff_analysis")
data(kostic2012crc)
set.seed(1024)
diffres <- diff_analysis(kostic2012crc, class="DIAGNOSIS",
		        mlfun="lda", filtermod="fdr",
		        firstcomfun = "kruskal.test", 
		        firstalpha=0.05, strictmod=TRUE,
	                submin=3, subclwilc=TRUE,
	                secondcomfun = "wilcox.test",
	                secondalpha=0.01, lda=3)

test_that("get_call",{
    expect_true(is(diffres@call, "call"))
})

context("get_call function")
test_that("get_call",{
    expect_equal("kruskal.test", MicrobiotaProcess:::get_call(diffres, "firstcomfun"))
    expect_equal("wilcox.test", MicrobiotaProcess:::get_call(diffres, "secondcomfun"))
    expect_equal("lda", MicrobiotaProcess:::get_call(diffres, "mlfun"))
})

