context("get_taxadf work")
test_that("get_taxadf work",{
    library(phyloseq)
    dir <- system.file("extdata", package="phyloseq")
    richph <- import_biom(file.path(dir,"rich_dense_otu_table.biom"))
    phyda <- get_taxadf(richph, taxlevel=2)
    expect_true(is(phyda, "phyloseq"))
})  
