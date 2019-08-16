context("gettaxdf work")
test_that("getaxdf work",{
    library(phyloseq)
    dir <- system.file("extdata", package="phyloseq")
    richph <- import_biom(file.path(dir,"rich_dense_otu_table.biom"))
    phyda <- gettaxdf(richph, taxlevel=2)
    expect_true(is(phyda, "phyloseq"))
})  
