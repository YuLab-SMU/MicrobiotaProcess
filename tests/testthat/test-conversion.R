context("convert_to_treedata")
test_that("data to treedata",{
    data(hmp_aerobiosis_small)
    treeda <- convert_to_treedata(taxda)
    expect_that(treeda, is_a("treedata"))
})
